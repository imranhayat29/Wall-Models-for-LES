%%**************************************************************************************%
%                         Integral wall model for LES                                   %
%***************************************************************************************%
% This code is based on the integral wall model formulation originally                  %
% developed and outlined by Yang et. al in the following paper:                         %
%                                                                                       %
% "Integral wall model for large eddy simulations of wall-bounded                       % 
%  turbulent flows." Physics of Fluids 27.2 (2015): 025112.                             %
%                                                                                       %
% Modifications to the original wall-model that you see in this code                    %
% are outlined in our upcoming work available on arXiv:                                 %
%                                                                                       %
%  https://arxiv.org/abs/2111.02542                                                     % 
%                                                                                       %
% "Implementation and performance analysis of efficient grid-free integral              % 
%  wall models in unstructured-grid LES solvers" (Hayat & Park, 2022)                   %
%                                                                                       %
%***************************************************************************************%
% Accompanying files required:                                                          %
%                                                                                       %   
% 1. Channel_DNS_JHTDB_Retau1000 (contains 2 files: mean_data_jhu.xlsx,re-tau.xlsx)     %                                                                                       %    
% 2. calc_grad.m                                                                        %
% 3. DNS_dt=0.0325_ndiv=16_yp100.mat (can be downloaded from the link below)            %
%                                                                                       %
%***************************************************************************************%
% Channel DNS data for the apriori testing of this code was acquired                    % 
% from Johns Hopkins Turbulence database. A sample DNS dataset is made                  %
% available on the following link:                                                      %
%                                                                                       %
% https://drive.google.com/drive/folders/1ESNdOrywOvNusEzFmAp7BWof_XUpvkUi?usp=sharing  %
%                                                                                       %
%***************************************************************************************%
% Questions/feedback can be directed to the following email:                            %               
% hayat3@seas.upenn.edu                                                                 %
%                                                                                       %
%***************************************************************************************%

% Note: x=streamwise , z=spanwise , y=wall-normal

% Variable naming from paper:
%
% Integral terms
%  Lx  : iwm_L(i,j,idx_Lu)
%  Lxx : iwm_L(i,j,idx_Luu)
%  Lz  : iwm_L(i,j,idx_Lw)
%  Lzz : iwm_L(i,j,idx_Lww)
%  Lxz : iwm_L(i,j,idx_Luw)
%   
% (tau_hwm-tau_w) : del_tau

clc
clear all

global idx_dirx idx_dirz idx_dim idx_Lu idx_Luu idx_Lw idx_Lww idx_Luw

%**************************************************************************
%% Define Indices
%**************************************************************************

% direction x (local streamwise)
idx_dirx = 1;
% direction z (local spanwise)
idx_dirz = 2;

% local dimensions of wall, wall model always deal with 2D surfaces 
idx_dim  = 2;

% integrated profile dimensions
idx_Lu  = 1;   % index for integral of u
idx_Luu = 2;   % index for integral of uu
idx_Lw  = 3;   % etc.
idx_Lww = 4;
idx_Luw = 5;

%**************************************************************************
%% Function INPUTS
%**************************************************************************

%========================================
% Constants:
%========================================
vonk  = 0.4;
B     = 5;

% Length of domain in streamwise and spanwise direction.
% (Needed for dx and dz for computing derivatives using Finite difference)
lx    = 8*pi;         
lz    = 3*pi;         

%========================================
% User defined:
%========================================
nt_end  = 30;          % Time step upto which you want to run the simulation 
                        % Last time step for accompanying JHUTDB DNS data 
                        % is 799
                        
theta   = 1;            % theta is the multiples of timescale used for time filtering 
iwm_tol = 0.000001;     % for Newton Raphson
iwm_eps = 0.000000001;  % for Newton Raphson
MaxIter = 1500;         % for Newton Raphson 

%Grid:
ndiv    = 16;           % for undersampled DNS data (specific to JHU test case)
nx      = 2048/ndiv;    % total x grid points for FVM etc...
nz      = 1536/ndiv;          
dx      = lx/nx;        % grid size in x-direction
dz      = lz/nz;        % grid size in spanwise-direction

% hwm     = 0.1008     
% dt      = 0.0065;     % hwm and LES timestep dt (for DNS test: database timestep)
                        % are commented out here because they are stored in
                        % the input DNS data file, based on at what height and
                        % timesteps the DNS data was extracted.

%=====================================
% DNS Data:
%=====================================                      
                        
%load DNS data
Ret_av_DNS = 999.35;    %fixed for JHTDB channel flow
ut_av_DNS  = 0.049968;  %fixed for JHTDB channel flow
nu         = 0.00005;   %fixed for JHTDB channel flow

% NOTE: The following .mat file contains JHTDB Channel DNS data from a plane at a
% wall-normal distance of yplus=100 or y=0.1008 from the bottom wall. The 
% data has been undersampled 16 times in each of the wall-parallel directions.
% The snapshots are collected every dt=0.0325. 

% The arrays contain data in the following form:
% [k,i,j] : k is time dimension, i is streamwise and j is spanwise.

load('DNS_dt=0.0325_ndiv=16_yp100.mat')  % Can be downloaded from the link 
                                         % given on the top in the description.
U_DNS   = U_full;
W_DNS   = W_full;
p_DNS   = P_full;
rho_DNS = ones(nt_end+1,nx,nz);          % DNS test case is incompressible flow
clear U_full W_full P_full;

% set the wall-model timestep as the same as LES timestep:
iwm_dt = dt;

% set h_wm for each wall face (here it is constant, since DNS data is taken
% from a plane at fixed height
iwm_hwm(1:nx,1:nz) = hwm;

%======================================
% Allocate size to arrays 
%======================================

% These variables are only Input to the WM:
u = zeros(nx,nz);
w = zeros(nx,nz);          
p = zeros(nx,nz);    
grad_p_filt = zeros(nx,nz,idx_dim);          % 2 components for each wall face
grad_L      = zeros(nx,nz,idx_Luw,idx_dim);  % 5 L's (Lu,Luu,Lw,Lww,Luw), 2
                                             % components for each.

% These variables are only output from the WM:
tauw_x = zeros(nx,nz,nt_end);
tauw_z = zeros(nx,nz,nt_end);

% These variables are both input/output to/from the WM:
U_tan_filt = zeros(nx,nz,idx_dim);
del_tau    = zeros(nx,nz,idx_dim);
iwm_L      = zeros(nx,nz,idx_Luw);
p_filt     = zeros(nx,nz);
iwm_filter = zeros(nx,nz);
utau_filt  = zeros(nx,nz);

%**************************************************************************
%% Function OUTPUTS
%**************************************************************************
 
% mean tau_w_x
% mean tau_w_z
% mean wallmodel velocity profiles

%**************************************************************************
%%                          MAIN ROUTINE                                  %
%**************************************************************************

%================================
% Initialization
%================================

% Initialize tau_w and all other quantities in the wall model
% (we use mean value of U-LES to approximate initial profiles at all grid
% points on wall)

Ui = mean(U_DNS(1,:,:),'all');       
Wi = mean(W_DNS(1,:,:),'all');
Pi = mean(p_DNS(1,:,:),'all');
rhoi = mean(rho_DNS(1,:,:),'all');

[tauw_x(:,:,1),tauw_z(:,:,1),p_filt,iwm_filter, ...
 utau_filt,del_tau,U_tan_filt,iwm_L] = initialize_iwm(nx,nz,hwm,vonk,B,nu,iwm_tol, ...
                                                  iwm_eps,Ui,Wi,Pi,rhoi,iwm_dt);

tic
%====================================================================
% Loops over all time steps available in the input field arrays
%====================================================================

for jt=1:nt_end
    
    fprintf('t=%5.4d (%5.2f%% complete)\n', jt*dt, 100*jt/nt_end)
    
    %Instantaneous DNS data (LES data in WMLES) from the matching plane (2 dimentional arrays)
    u(:,:)  = U_DNS(jt+1,:,:);
    w(:,:)  = W_DNS(jt+1,:,:);   
    p(:,:)  = p_DNS(jt+1,:,:);
    rho(:,:)= rho_DNS(jt+1,:,:);
    
    % Compute the surface gradients required on RHS of eqn (C22) in Yang et al 2015
    % Note that each gradient has 2 components, one in x and the other in z.
    % The subroutine calc_grad computes the surface gradients with periodic
    % boundary conditions assumed.
    
    [grad_L(:,:,idx_Lu ,idx_dirx), grad_L(:,:,idx_Lu ,idx_dirz) ]= calc_grad(nx,nz,dx,dz,iwm_L(:,:,idx_Lu) );
    [grad_L(:,:,idx_Luu,idx_dirx), grad_L(:,:,idx_Luu,idx_dirz) ]= calc_grad(nx,nz,dx,dz,iwm_L(:,:,idx_Luu));
    [grad_L(:,:,idx_Lw ,idx_dirx), grad_L(:,:,idx_Lw ,idx_dirz) ]= calc_grad(nx,nz,dx,dz,iwm_L(:,:,idx_Lw) );
    [grad_L(:,:,idx_Lww,idx_dirx), grad_L(:,:,idx_Lww,idx_dirz) ]= calc_grad(nx,nz,dx,dz,iwm_L(:,:,idx_Lww));    
    [grad_L(:,:,idx_Luw,idx_dirx), grad_L(:,:,idx_Luw,idx_dirz) ]= calc_grad(nx,nz,dx,dz,iwm_L(:,:,idx_Luw));    
    [grad_p_filt(:,:,idx_dirx)   , grad_p_filt(:,:,idx_dirz)    ]= calc_grad(nx,nz,dx,dz,p_filt(:,:));

    
    % NOTE: At this point all required inputs to the wall model are 
    % available at each wall face.    
    
    %======================================================================
    % Loop over all wall faces and compute the wall stress at each wall face
    %======================================================================      
    for iwm_i = 1:nx
        for iwm_j = 1:nz
            
            % local values at each wall face to be fed to iwm_calc_wallstress routine
            hwm_l        = iwm_hwm     (iwm_i,iwm_j);
            rho_l        = rho         (iwm_i,iwm_j);
            u_l          = u           (iwm_i,iwm_j);
            w_l          = w           (iwm_i,iwm_j);
            p_l          = p           (iwm_i,iwm_j);
            p_filt_l     = p_filt      (iwm_i,iwm_j);
            iwm_filter_l = iwm_filter  (iwm_i,iwm_j);
            utau_filt_l  = utau_filt   (iwm_i,iwm_j);           
            grad_p_l     = grad_p_filt (iwm_i,iwm_j,:);
            del_tau_l    = del_tau     (iwm_i,iwm_j,:);
            U_tan_filt_l = U_tan_filt  (iwm_i,iwm_j,:);
            iwm_L_l      = iwm_L       (iwm_i,iwm_j,:);
            grad_L_l(:,:)= grad_L      (iwm_i,iwm_j,:,:);
            
            %==============================================================
            % Call the wall-stress calculation routine to get wall stress at
            % the current wallface at (iwm_i,iwm_j)
            %==============================================================             
            
            % the subroutine to calculate wall stress
            [twx,twz,Ax,Cx,p_filt_l,iwm_filter_l,...
             utau_filt_l,del_tau_l,U_tan_filt_l, ...
             iwm_L_l,equil_flag] ...
                     = iwm_calc_wallstress(hwm_l,rho_l,u_l,w_l,p_l,p_filt_l, ...
                       iwm_tol,iwm_eps,MaxIter,nu,vonk,B,theta,grad_L_l,grad_p_l,  ...
                       iwm_i,iwm_j,iwm_dt,iwm_filter_l,utau_filt_l,del_tau_l,...
                       U_tan_filt_l,iwm_L_l);

            % Imposing tauw_x, tauw_z in the LES solver as a boundary
            % condition wall stress
            tauw_x(iwm_i,iwm_j,jt) = twx;
            tauw_z(iwm_i,iwm_j,jt) = twz;

            % Update variables that will be used by the WM at next timestep                   
            p_filt    (iwm_i,iwm_j)   =     p_filt_l;
            iwm_filter(iwm_i,iwm_j)   = iwm_filter_l;            
            utau_filt (iwm_i,iwm_j)   =  utau_filt_l;
            del_tau   (iwm_i,iwm_j,:) =    del_tau_l;
            U_tan_filt(iwm_i,iwm_j,:) = U_tan_filt_l;
            iwm_L     (iwm_i,iwm_j,:) =      iwm_L_l;
            
            %==============================================================
            % Store velocity profiles for plotting
            %==============================================================              
            
            % Instantaneous velocity profile
            utx = sign(twx)*sqrt(abs(twx)/rho_l);
            utz = sign(twz)*sqrt(abs(twz)/rho_l);
            
            Vel = sqrt(U_tan_filt_l(idx_dirx)^2. +U_tan_filt_l(idx_dirz)^2.);
            utau = (utx^4.+utz^4.)^0.25;

            yy      = logspace(log10(11*nu/utau),log10(1),100);
            yy_shrt = logspace(log10(0.0001), log10(11*nu/utau), 50);
            
            if equil_flag == 0
                del_visc = nu*sqrt(utx^2.+utz^2.)/(utau^2.0);
                
                ulog(iwm_i,iwm_j,:)= utau*((U_tan_filt_l(idx_dirx)/Vel)*...
                                (log(yy/hwm_l)/vonk)  + Ax*(yy/hwm_l) + Cx);
                uvisc(iwm_i,iwm_j,:)=utx*yy_shrt/del_visc;
            else
                del_visc = nu/utau;
                
                ulog(iwm_i,iwm_j,:)= utau*(U_tan_filt_l(idx_dirx)/Vel)*...
                                          (log(yy/del_visc)/vonk + B);
                uvisc(iwm_i,iwm_j,:)=utau*(U_tan_filt_l(idx_dirx)/Vel)*...
                                          yy_shrt/del_visc;
            end           
        end
    end
    
    % Collect some Stats
    tauw_tot(:,:,jt) = sqrt(tauw_x(:,:,jt).^2 + tauw_z(:,:,jt).^2);
    tauw_wm(jt)      = mean(tauw_tot(:,:,jt),'all');

    % Spatially averaged velocity profile
    ulog_m(jt,:)  = mean(ulog ,[1 2]);
    uvisc_m(jt,:) = mean(uvisc,[1 2]);   
    
end
toc

% Mean (time-averaged) velocity profiles
ulog_p = mean(ulog_m(10:end,:),1)/ut_av_DNS;
uvisc_p = mean(uvisc_m(10:end,:),1)/ut_av_DNS;

%**************************************************************************
%% PLOT RESULTS
%**************************************************************************
 
%===============================
% Plot Profiles
%===============================

%Import DNS data for plotting profile
ll     = xlsread('Channel_DNS_JHTDB_Retau1000/mean_data_jhu');   
ypll   = ll(:,2);
Upll   = ll(:,3);
semilogx(ypll,Upll,'k.','MarkerSize',5);
hold on

semilogx(yy(1:49)/(nu/ut_av_DNS),ulog_p(1:49),'r','LineWidth',2.5);
semilogx(yy_shrt/(nu/ut_av_DNS),uvisc_p,'b','LineWidth',2.5);

xlabel('$y^{+}$','Interpreter','Latex');
ylabel('$U^{+}$','Interpreter','Latex');
xlim([0.1 2000]);
ylim([0 26]);


legend('mean-DNS (JHU Re_{\tau}=1000)','log-layer','viscous sublayer','Location','NorthWest')

%================================
% Plot tau_w ratio vs time
%================================
mm       = xlsread('Channel_DNS_JHTDB_Retau1000/re-tau');
time_dns = mm(:,1);
ret_dns  = mm(:,2);
tw_DNS   = (nu*ret_dns).^2;

tt = 1:nt_end;
tratio_wm  = abs(tauw_wm)/ut_av_DNS^2;
tratio_DNS = tw_DNS/ut_av_DNS^2;

figure
plot(tt*dt,tratio_wm,'r');
hold on
plot(time_dns,tratio_DNS,'k');
xlabel('t');
xlim([0 26])
ylim([0.95 1.05]);
ylabel('\tau_w/\tau_{w,DNS av}')

legend('IWM','DNS','Location','NorthEast')
title('Spatially-averaged wall stress magnitude normalized by mean DNS wall stress')



%**************************************************************************
%%                           SUBROUTINES                                   %

%**************************************************************************
%% Subroutine: initialize_iwm
%**************************************************************************
% This subroutine initializes the wall model with plug flow conditions

function [tauw_x,tauw_z,p_filt,iwm_filter,utau_filt,del_tau, ...
         U_tan_filt,iwm_L] = initialize_iwm(nx,nz,hwm,vonk,B,nu,iwm_tol, ...
                                        iwm_eps,Ui,Wi,Pi,rhoi,dt)

global idx_dirx idx_dirz idx_dim idx_Lu idx_Luu idx_Lw idx_Lww idx_Luw
        
% initial value for the x and z-velocity at first grid point
uinit = Ui;
winit = Wi;

% filitered velocity at the first grid point in x, z directions
U_tan_filt(1:nx,1:nz,idx_dirx) = uinit;
U_tan_filt(1:nx,1:nz,idx_dirz) = winit;

%==========================================================================
% Note: We initialize utx and utz with standard log-law for equilibrium 

Vel = sqrt(uinit^2. + winit^2.);
Uproj = uinit/Vel;
Wproj = winit/Vel;

equilut = sign(Vel);
feqm = (equilut/vonk)*log(hwm*equilut/nu) + B*equilut - Vel;

iterm=0;
while abs(feqm) > iwm_tol && iterm <= 20
    ut_p = equilut + iwm_eps;
    feqm_ps = (ut_p/vonk)*log(hwm*ut_p/nu) + B*ut_p - Vel;
    a = (feqm_ps - feqm)/iwm_eps;
    equilut = equilut - feqm/a;
    
    feqm = (equilut/vonk)*log(hwm*equilut/nu) + B*equilut - Vel;
    iterm = iterm + 1;
end

equilutx = real(equilut*sqrt(Uproj));
equilutz = real(equilut*sqrt(Wproj));

%==========================================================================

utau = (equilutx^4.+equilutz^4.)^0.25;
deli = 11*nu/utau;
delv = nu/utau;
rat  = deli/hwm;

% wall stress in x, z directions
tauw_x(1:nx,1:nz) = sign(equilutx)*rhoi*equilutx^2.;
tauw_z(1:nx,1:nz) = sign(equilutz)*rhoi*equilutz^2.;

% pressure at first grid point
p_filt(1:nx,1:nz) = Pi;

% integrals of Lu, Lw, etc.
iwm_L(1:nx,1:nz,idx_Lu)  = utau*Uproj*(0.5*(deli^2.0)/delv + hwm*( B*(1-rat)...
                              - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(hwm/delv)) ));
        
iwm_L(1:nx,1:nz,idx_Lw)  = utau*Wproj*(0.5*(deli^2.0)/delv + hwm*( B*(1-rat)...
                              - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(hwm/delv)) ));
                          
iwm_L(1:nx,1:nz,idx_Luu) = (utau*Uproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
            hwm*((1.0/vonk)^2*((log(hwm/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
            +(2*B/vonk)*((log(hwm/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat))); 
        
iwm_L(1:nx,1:nz,idx_Lww) = (utau*Wproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
            hwm*((1.0/vonk)^2*((log(hwm/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
            +(2*B/vonk)*((log(hwm/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));  
        
iwm_L(1:nx,1:nz,idx_Luw) = (utau^2)*Uproj*Wproj*(1.0/3.0*(deli^3.0/delv^2.0)+...
            hwm*((1.0/vonk^2)*((log(hwm/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
            +(2*B/vonk)*((log(hwm/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));  
   
% RHS terms in the integral equation to be used in the next timestep
del_tau(1:nx,1:nz,1:idx_dim) = 0.;

% filtered friction velocity and the filtering time scale
utau_filt(1:nx,1:nz) = utau;
iwm_filter(1:nx,1:nz)= dt/(hwm/vonk/equilutx); 

end

%**************************************************************************
%% Subroutine: iwm_calc_wallstress
%**************************************************************************

function [tauw_x,tauw_z,Ax,Cx,p_filt,iwm_filter,utau_filt,del_tau,U_tan_filt, ...
          iwm_L,equil_flag] = iwm_calc_wallstress(hwm,rho,u,v,p,p_filt, ...
          iwm_tol,iwm_eps,MaxIter,nu,vonk,B,theta, grad_L,grad_p_filt, ...
          iwm_i,iwm_j,iwm_dt,iwm_filter,utau_filt,del_tau,U_tan_filt,iwm_L)

global idx_dirx idx_dirz idx_Lu idx_Luu idx_Lw idx_Lww idx_Luw
         
% Externally imposed mean streamwise pressure gradient.
dpdx_imposed = 0;   

%==========================================================================
% Calculation of RHS at timestep (n-1)
%==========================================================================
% This section calculates the right hand side of the iwm system in 
% eqn (C22) in Yang et al 2015.

% Note that this is named as iwm_lhs because all the terms on the RHS of 
% eqn (C22) are taken to the LHS, thus the reversed sign for each term 
% compared to eqn (C22).

% Note that in calculating RHS, all the terms are taken from the previous
% timestep (n-1).

    % convective terms	
    conv_x = grad_L(idx_Luu,idx_dirx) + grad_L(idx_Luw,idx_dirz) - ...
             U_tan_filt(idx_dirx) * ( grad_L(idx_Lu,idx_dirx)  + ...
                                      grad_L(idx_Lw,idx_dirz) );
                                  
    conv_z = grad_L(idx_Luw,idx_dirx) + grad_L(idx_Lww,idx_dirz) - ...
             U_tan_filt(idx_dirz) * ( grad_L(idx_Lu,idx_dirx)  + ...
                                      grad_L(idx_Lw,idx_dirz) );			

    % Include the mean pressure gradient if any
    grad_p_filt(idx_dirx) =  grad_p_filt(idx_dirx) - dpdx_imposed;	

    % RHS term:
    % This is LHS in discretized ODE: 
    % y^(n)-y^(n-1)-dt*f(y,t)^(n-1) = y^(n) + LHS = 0, where y = Lx,Lxx etc. 
    % OR the integrated momentum equation (C22), except for the Lu term
    iwm_lhs_x = -iwm_L(idx_Lu) + iwm_dt*( conv_x + (1.0/rho)* ...
                (grad_p_filt(idx_dirx) * hwm - del_tau(idx_dirx)) );	
    iwm_lhs_z = -iwm_L(idx_Lw) + iwm_dt*( conv_z + (1.0/rho)* ...
                (grad_p_filt(idx_dirz) * hwm - del_tau(idx_dirz)) ); 
          
%==========================================================================
% temporal filtering at (n)
%==========================================================================
    U_tan_filt(idx_dirx) = U_tan_filt(idx_dirx)*(1.0 - iwm_filter) + u * iwm_filter;
    U_tan_filt(idx_dirz) = U_tan_filt(idx_dirz)*(1.0 - iwm_filter) + v * iwm_filter;
    p_filt = p_filt*(1.0 - iwm_filter) + p * iwm_filter;
    
% Note that all variables after this point are at current timestep (n).
 
%==========================================================================
% Equilibrium WM solution 
%(for initialization of utau and reverting back, in case of failure in iWM)
%========================================================================== 

    % Define some useful flags:
    deli_flag  = 0;            % activates if hwm^+ <= 11.
    equil_flag = 0;            % activates when iWM solution fails to converge.

% Note: We initialize utx and utz with standard log-law used in equilibrium flag 

    %total velocity
    Vel = sqrt(U_tan_filt(idx_dirx)^2. + U_tan_filt(idx_dirz)^2.);
    
    Uproj = U_tan_filt(idx_dirx)/Vel;
    Wproj = U_tan_filt(idx_dirz)/Vel;

    equilut = sign(Vel);
    
    feqm = (equilut/vonk)*log(hwm*equilut/nu) + B*equilut - Vel;
    
    % use Newton method to solve the log law.
    iterm=0;
    while abs(feqm) > iwm_tol && iterm <= 20
        a = (1.0/vonk)*(log(hwm*equilut/nu) + 1.0) + B;        
        equilut = equilut - feqm/a;
        
        feqm = (equilut/vonk)*log(hwm*equilut/nu) + B*equilut - Vel;
        iterm = iterm + 1;
    end   
    
    if iterm==21
        fprintf('ut_eqm did not converge in 20 iterations for i=%i j=%i\n',iwm_i,iwm_j);
    end 
    
    
    if abs(hwm) <= abs(11*nu/equilut)
        equilut = sqrt(nu*Vel/hwm);
        deli_flag=1;
    end
    
    equilutx = sign(Uproj)*equilut*sqrt(abs(Uproj));
    equilutz = sign(Wproj)*equilut*sqrt(abs(Wproj));

%==========================================================================
% Integral WM solution
%==========================================================================    
    
    %Seed values of utx & utz
    utau_x = equilutx;
    utau_z = equilutz;
    
    % iwm_slv evaluates eqns (C23) and (C24) for the given values of utx, 
    % utz and RHS (iwm_lhs). The values of fx and fz should be close to
    % zero (i.e. < iwm_tol) for the converged values of utx and utz.
        
    [fx,fz] = iwm_slv(iwm_lhs_x,iwm_lhs_z,U_tan_filt(idx_dirx), ...
                      U_tan_filt(idx_dirz),hwm, utau_x, utau_z, nu, vonk);   

    iter = 0;   
    
    if deli_flag==0
        
        % use Newton method to solve the system
        while max(abs(fx),abs(fz)) > iwm_tol
            
            iwmutxP = utau_x+iwm_eps;
            iwmutzP = utau_z;
            
            [fxp, fzp] = iwm_slv(iwm_lhs_x,iwm_lhs_z,U_tan_filt(idx_dirx),...
                         U_tan_filt(idx_dirz), hwm, iwmutxP, iwmutzP, nu, vonk);
            
            a11 = (fxp-fx)/iwm_eps;
            a21 = (fzp-fz)/iwm_eps;
            
            iwmutxP = utau_x;
            iwmutzP = utau_z+iwm_eps;
            
            [fxp, fzp] = iwm_slv(iwm_lhs_x,iwm_lhs_z,U_tan_filt(idx_dirx),...
                         U_tan_filt(idx_dirz), hwm, iwmutxP, iwmutzP, nu, vonk);
            
            a12 = (fxp-fx)/iwm_eps;
            a22 = (fzp-fz)/iwm_eps;
            
            utau_x = utau_x - 0.50*( a22*fx-a12*fz)/(a11*a22-a12*a21);
            utau_z = utau_z - 0.50*(-a21*fx+a11*fz)/(a11*a22-a12*a21);
            
            % infinity check
            if (abs(utau_x) > 10^4.0 || abs(utau_z)> 10^4.0)
                fprintf('divergence in utx or utz at i=%i j=%i , reverting to Equilibrium WM \n',iwm_i,iwm_j);
                equil_flag = 1;
                break
            end
            
            [fx, fz] = iwm_slv(iwm_lhs_x, iwm_lhs_z,U_tan_filt(idx_dirx), ...
                            U_tan_filt(idx_dirz),hwm, utau_x, utau_z, nu, vonk);
            
            iter = iter+1;
            
            % maximum iteration reached
            if (iter>MaxIter)
                equil_flag = 1;
                break
            end
        end
    end

%==========================================================================
% Switch to Equilibrium WM if eqm_flag is 1
%==========================================================================  
 
    if (equil_flag==1 || deli_flag==1)
        utau_x = equilutx;
        utau_z = equilutz;
    end

%==========================================================================
% Evaluate all other quantities (based on the converged values of utx & utz)
%==========================================================================    
    
    % At this point utx and utz have been calculated. After this step, we simply
    % calculate other parameters using their relations with utx, utz and utau
    
    % calculate the friciton velocity
    utau = (utau_x^4.+utau_z^4.)^0.25;
        
    % eq. C25 in Yang et al. 2015
    deli = min( 11*nu/utau , hwm ); 
    
    % It is convenient to define this ratio to avoid clutter in eqns below.
    rat = deli/hwm;
    
    %calculate Ax, Az, Cx, Cz
    if(equil_flag==1  || deli_flag==1)
        
        %By definition of single profile characterized by utau for eqm flag
        delv = nu/utau;
        
        Ax = 0.;
        Az = 0.;
        Cx = 0.;
        Cz = 0.;
    else
        % Reduced form of eqn. (C26) in Yang et al. 2015 when modification 
        % proposed in Hayat & Park 2021 is applied to the viscous profile 
        % (see eqns (2.30) and (A5) in Hayat & Park 2021)
        delv = nu/utau;
        
        if deli == 11*nu/utau
            % eq. C27 in Yang et al. 2015 with Hayat & Park 2021 modification
            Ax = (U_tan_filt(idx_dirx)/utau + Uproj/vonk*log(rat) - ...
                (deli/delv)*sign(utau_x)*(utau_x/utau)^2.0 ) / ((1.0-rat));
            Az = (U_tan_filt(idx_dirz)/utau + Wproj/vonk*log(rat) - ...
                (deli/delv)*sign(utau_z)*(utau_z/utau)^2.0 ) / ((1.0-rat));
            Cx = U_tan_filt(idx_dirx)/utau - Ax;
            Cz = U_tan_filt(idx_dirz)/utau - Az;
            
        elseif deli == hwm
            Ax = 0.;
            Az = 0.;
            Cx = 0.;
            Cz = 0.;
            
        end
    end

    % check for excessive linear term correction
    % Revert to EQWM if values exceed 2.0
    if (abs(Ax)>2|| abs(Az)>2)
        equil_flag = 1;
        
        utau_x = equilutx;
        utau_z = equilutz;
        Ax = 0.;
        Az = 0.;
        Cx = 0.;
        Cz = 0.;
        
        utau = (utau_x^4.+utau_z^4. )^0.25;
        deli = 11*nu/utau;
        
        if deli_flag==1
            deli = hwm;
        end
        
        delv = nu/utau;
        rat  = deli/hwm;        
    end
   
    % compute the required integrals   
    
    if (equil_flag==1 || deli_flag==1)
        
        % Equilibrium solution:
        
        if deli_flag==0
            % Lu
            iwm_L(idx_Lu) = utau*Uproj*(0.5*(deli^2.0)/delv + ...
                hwm*( B*(1-rat) - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(hwm/delv)) ));
            % Lw
            iwm_L(idx_Lw) = utau*Wproj*(0.5*(deli^2.0)/delv + ...
                hwm*( B*(1-rat) - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(hwm/delv)) ));
            % Luu
            iwm_L(idx_Luu) = (utau*Uproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
                hwm*((1.0/vonk)^2*((log(hwm/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
                +(2*B/vonk)*((log(hwm/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));
            % Lww
            iwm_L(idx_Lww) = (utau*Wproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
                hwm*((1.0/vonk)^2*((log(hwm/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
                +(2*B/vonk)*((log(hwm/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));
            % Luw
            iwm_L(idx_Luw) = (utau^2)*Uproj*Wproj*(1.0/3.0*(deli^3.0/delv^2.0)+...
                hwm*((1.0/vonk^2)*((log(hwm/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
                +(2*B/vonk)*((log(hwm/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));
        else
            % Lu:  Eq. C19 in Yang et al. 2015 with Hayat & Park 2021 modification
            iwm_L(idx_Lu) = 0.5*sign(utau_x)*((utau_x^2.0)/utau)*(deli^2.0)/delv;
            % Lw
            iwm_L(idx_Lw) = 0.5*sign(utau_z)*((utau_z^2.0)/utau)*(deli^2.0)/delv;
            
            % Luu: Eq. C20 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_L(idx_Luu) = 1.0/3.0*(utau_x^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            % Lww
            iwm_L(idx_Lww) = 1.0/3.0*(utau_z^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            
            % Luw: Eq. C21 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_L(idx_Luw) = 1.0/3.0*sign(utau_x)*sign(utau_z)*(utau_x*utau_z/utau)^2.0*(deli^3.0/delv^2.0);
        end

    else
        % Integral WM solution:
        
        if deli == 11*nu/utau
            
            % Lu:  Eq. C19 in Yang et al. 2015 with Hayat & Park 2021 modification
            iwm_L(idx_Lu) = 0.5*sign(utau_x)*((utau_x^2.0)/utau)*(deli^2.0)/delv + ...
                utau*hwm*(0.5*Ax*(1-rat^2.0) + Cx*(1-rat) - Uproj/vonk*(1-rat+rat*log(rat)));
            % Lw
            iwm_L(idx_Lw) = 0.5*sign(utau_z)*((utau_z^2.0)/utau)*(deli^2.0)/delv + ...
                utau*hwm*(0.5*Az*(1-rat^2.0) + Cz*(1-rat) - Wproj/vonk*(1-rat+rat*log(rat)));
            
            % Luu: Eq. C20 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_L(idx_Luu) = 1.0/3.0*(utau_x^4.0/utau^2.0)*(deli^3.0/delv^2.0)+...
                utau^2.0*hwm*(-Ax*(Uproj/vonk)*rat^2.0*log(rat) + Ax*(Cx-Uproj/vonk/2.0) ...
                *(1-rat^2.0) + Ax^2.0/3.0*(1-rat^3.0) + (Cx - Uproj/vonk)^2.0   ...
                -rat*(Cx-Uproj/vonk+Uproj/vonk*log(rat))^2.0 + (Uproj/vonk)^2.0*(1-rat));
            % Lww
            iwm_L(idx_Lww) = 1.0/3.0*(utau_z^4.0/utau^2.0)*(deli^3.0/delv^2.0)+...
                utau^2.0*hwm*(-Az*(Wproj/vonk)*rat^2.0*log(rat) + Az*(Cz-Wproj/vonk/2.0) ...
                *(1-rat^2.0) + Az^2.0/3.0*(1-rat^3.0) + (Cz - Wproj/vonk)^2.0   ...
                -rat*(Cz-Wproj/vonk+Wproj/vonk*log(rat))^2.0 + (Wproj/vonk)^2.0*(1-rat));
            
            % Luw: Eq. C21 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_L(idx_Luw) = 1.0/3.0*sign(utau_x)*sign(utau_z)*(utau_x*utau_z/utau)^2.0*(deli^3.0/delv^2.0)+...
                utau^2.0*hwm*(-1.0/vonk*(Ax*Wproj+Az*Uproj)*(0.25-0.25*rat^2.0+0.5*rat^2.0*log(rat))...
                -1/vonk*(Cx*Wproj+Cz*Uproj)*(1-rat+rat*log(rat))-Uproj*Wproj/(vonk^2.0)*(rat-2+rat*(log(rat)-1)^2.0)...
                +(1.0/3.0)*Ax*Az*(1-rat^3.0)+0.5*(Ax*Cz+Az*Cx)*(1-rat^2.0)+Cx*Cz*(1-rat));
            
        elseif deli == hwm
            
            % Lu:  
            iwm_L(idx_Lu) = 0.5*sign(utau_x)*((utau_x^2.0)/utau)*(deli^2.0)/delv;
            % Lw
            iwm_L(idx_Lw) = 0.5*sign(utau_z)*((utau_z^2.0)/utau)*(deli^2.0)/delv;            
            % Luu:
            iwm_L(idx_Luu) = 1.0/3.0*(utau_x^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            % Lww
            iwm_L(idx_Lww) = 1.0/3.0*(utau_z^4.0/utau^2.0)*(deli^3.0/delv^2.0);            
            % Luw: 
            iwm_L(idx_Luw) = 1.0/3.0*sign(utau_x)*sign(utau_z)*(utau_x*utau_z/utau)^2.0*(deli^3.0/delv^2.0);
        end
    end

    if  deli_flag==1
        
        % calculate top derivatives
        dudyT_x = U_tan_filt(idx_dirx)/hwm;
        dudyT_z = U_tan_filt(idx_dirz)/hwm;
        
        % calculte the turbulent diffusion term
        dVeldyT = Uproj*U_tan_filt(idx_dirx)/hwm + Wproj*U_tan_filt(idx_dirz)/hwm;        

    else
        
        % calculate top derivatives
        % From Eq. C12 Yang et al. 2015
        dudyT_x = (Uproj/vonk+Ax)*utau/hwm;
        dudyT_z = (Wproj/vonk+Az)*utau/hwm;
        
        % calculte the turbulent diffusion term
        % Eq. C14 Yang et al. 2015
        % Note: this reduces to utau/(vonk*Dz) for eqm flag
        dVeldyT = (Uproj*(Uproj/vonk+Ax) + Wproj*(Wproj/vonk+Az))*utau/hwm;        
    end
               
    % Eq. C13 Yang et al. 2015, the eddy viscosity 'nut' at y=hwm is:
    nut_hwm = (vonk*hwm)^2.0*dVeldyT;
    
    % calculate the wall stress
    %Note: Unlike prescribed roughness case, here we have the same form for
    %tauwx for both equilibrium and non-equilibrium
    
    tauw_x = sign(utau_x)*rho*utau_x^2.0;
    tauw_z = sign(utau_z)*rho*utau_z^2.0;
    
    % shear stress difference between top and bottom
    del_tau(idx_dirx) = (nut_hwm+nu)*dudyT_x - tauw_x;   
    del_tau(idx_dirz) = (nut_hwm+nu)*dudyT_z - tauw_z;
    
    % TIMESCALE: the filtered friction velocity used for filtering time scale
    utau_filt = iwm_filter*utau + (1.0-iwm_filter)*utau_filt;
    
    % update the filtering time scale (epsilon in Eq. 26 Yang et al. 2015)    
    iwm_filter = iwm_dt/(theta*hwm/utau_filt/vonk);
    
    % filtering time scale can only be larger than the time step,
    % if not, then just use the instantaneous flow field to do the model
    if (iwm_filter > 1.)
        iwm_filter = 1.;
    end

end
%*******************************************************************************
%% Subroutine: iwm_slv
%*******************************************************************************
% iwm_slv evaluates eqns (C23) and (C24) in Yang et al. 2015 for the given 
% values of utx, utz and RHS (iwm_lhs). The values of fx and fz should be 
% close to zero (i.e. < iwm_tol) for the converged values of utx and utz.
    
function [fx,fz]= iwm_slv(lhsx,lhsz,Ux,Uz,hwm,utx,utz,nu,vonk)

%U_tot
Vel = sqrt(Ux^2.0+Uz^2.0);

%Eqn C30
utau = (utx^4.0 + utz^4.0)^(0.25);

%Eqn C26 (viscous scale)
delv = nu/utau;

%Eqn C25 
%(Interface of log and viscous layer chosen as min of y+=11 or h_wm)
deli = min(11*nu/utau,hwm);
% deli = 11*nu/utau;             % 11*delv; 

%Find parameters of log profile only if deli=11*nu/utau, otherwise linear
%profile in whole inner layer
if deli == 11*nu/utau
    %Eqn C27
    Ax = (Ux/utau + (Ux/Vel)/vonk*log(deli/hwm) - (deli/delv)*sign(utx)* ...
         (utx/utau)^2.0) / ((1.0-deli/hwm));
    Az = (Uz/utau + (Uz/Vel)/vonk*log(deli/hwm) - (deli/delv)*sign(utz)* ...
         (utz/utau)^2.0) / ((1.0-deli/hwm));
    
    %Eqn C28
    Cx = Ux/utau - Ax;
    Cz = Uz/utau - Az;
        
    %define ratio of deli to h_wm
    rat = (deli/hwm);
    
    %LHS of Eqn C23 and C24
    Lu = 0.5*sign(utx)*(utx^2.0/utau)*(deli^2.0)/delv + utau*hwm* ...
        ( 0.5*Ax*(1-rat^2.0) + Cx*(1-rat) - (Ux/Vel)/vonk*(1-rat+rat*log(rat)) );
    Lw = 0.5*sign(utz)*(utz^2.0/utau)*(deli^2.0)/delv + utau*hwm* ...
        ( 0.5*Az*(1-rat^2.0) + Cz*(1-rat) - (Uz/Vel)/vonk*(1-rat+rat*log(rat)) );
    
elseif deli == hwm
    
    Lu = 0.5*sign(utx)*(utx^2.0/utau)*(deli^2.0)/delv;
    Lw = 0.5*sign(utz)*(utz^2.0/utau)*(deli^2.0)/delv; 
    
end

fx = Lu+lhsx;
fz = Lw+lhsz; 

end