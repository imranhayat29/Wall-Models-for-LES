%%**************************************************************************************%
%                         Integral wall model for LES                                   %
%***************************************************************************************%
% This code is based on the integral wall model formulation developed and               %
% outlined originally by Yang et. al in the following paper:                            %
%                                                                                       %
% "Integral wall model for large eddy simulations of wall-bounded                       % 
%  turbulent flows." Physics of Fluids 27.2 (2015): 025112.                             %
%                                                                                       %
% Modifications to the original wall-model that you see in this code                    %
% are outlined in our upcoming work available on arXiv:                                 %
%                                                                                       %
% "Numerical implementation of efficient grid-free integral wall models                 %
%  in unstructured-grid LES solvers" (Hayat & Park, 2021)                               %
%                                                                                       %
%***************************************************************************************%
% DNS data for apriori testing of this code can be found at the following               %
% link:                                                                                 %
%                                                                                       %
% https://drive.google.com/drive/folders/1ESNdOrywOvNusEzFmAp7BWof_XUpvkUi?usp=sharing  %
%                                                                                       %
%***************************************************************************************%

% Note: x=streamwise , y=spanwise , z=wall-normal

clc
clear all

global iwm_dirx iwm_diry iwm_DN iwm_Lu iwm_Luu iwm_Lv iwm_Lvv iwm_Luv ...
       iwm_utx iwm_uty iwm_tauwx iwm_tauwy iwm_flt_tagvel iwm_flt_tagvel_m iwm_flt_p  ...
       iwm_conv iwm_PrsGrad iwm_diff iwm_lhs iwm_flt_us ...
       iwm_tR iwm_Dz iwm_Ax iwm_Ay iwm_Cx iwm_Cy iwm_dt iwm_inte iwm_inte_m ...
       equil_flag div_flag neg_flag A_flag maxiter_flag neg_flag_uty neg_flag_V

%**************************************************************************
%% Define Indices
%**************************************************************************

% direction x (streamwise)
iwm_dirx = 1;
% direction y (spanwise)
iwm_diry = 2;

% local dimensions of wall, wall model always deal with 2D surfaces 
iwm_DN  = 2;

% integrated profile dimensions
iwm_Lu  = 1;   % index for integral of u
iwm_Luu = 2;   % index for integral of uu
iwm_Lv  = 3;   % etc.
iwm_Lvv = 4;
iwm_Luv = 5;

%**************************************************************************
%% Function INPUTS
%**************************************************************************
%Constants:
vonk    = 0.4;
B       = 5;
L_x     = 8*pi;         % Length of domain in axial direction!
                        % Needed for dx for computing derivatives using Finite
                        % difference
L_y     = 3*pi;
L_z     = 1; %2;

%User defined
nt_end  = 50;          % Time step upto which you want to run the 
                        % simulation (Last time step for accompanying
                        % JHUTDB DNS data is 800

theta   = 1;            % theta is the multiples of timescale used for time filtering 
iwm_tol = 0.000001;     % for Newton Raphson
iwm_eps = 0.000000001;  % for Newton Raphson
MaxIter = 1500;         % for Newton Raphson 

%Grid:
ndiv    = 16;           % for undersampled DNS data (specific to JHU test case)
nx      = 2048/ndiv;    % total x grid points for FVM etc...
ny      = 1536/ndiv;          
nz      = 1; %512;
dx      = L_x/nx;       % grid size in x-direction
dy      = L_y/ny;       % grid size in spanwise-direction

% hwm     = 0.1008     
% dt      = 0.0065;     % hwm and LES timestep dt (for DNS test: database timestep)
                        % are commented out here because I had stored it in
                        % the input data file based on at what height and
                        % timesteps the DNS data was extracted.
                        
cfl     = 1;            % CFL of LES (ONLY USED IN IC)

%load DNS data
Ret_av_DNS = 999.35;    %fixed for JHTDB channel flow
ut_av_DNS  = 0.049968;  %fixed for JHTDB channel flow
nu         = 0.00005;   %fixed for JHTDB channel flow

% The following .mat file contains JHTDB Channel DNS data from a plane at a
% wall-normal distance of yplus=100 or y=0.1008 from the bottom wall. The 
% data has been undersampled 16 times in each of the wall-parallel directions.
% The snapshots are collected every dt=0.0325. 

% The arrays contain data in the following form:
% [k,i,j] : k is time dimension, i is streamwise and j is spanwise.

load('DNS_dt=0.0325_ndiv=16_yp100.mat')  % Can be downloaded from the link on the top.
U_DNS   = U_full;
V_DNS   = W_full;
W_DNS   = V_full;
p_DNS   = P_full;
clear U_full V_full W_full P_full;

jt_end  = nt_end;              % LES counter end-time
% jt_end  = size(U_DNS,1);       % LES counter end-time

u = zeros(nx,ny);
v = zeros(nx,ny);          
w = zeros(nx,ny);
p = zeros(nx,ny);

dz = 2*hwm;
%**************************************************************************
%% Function OUTPUTS
%**************************************************************************
%mean:
%   taux
%   tauy
%   wallmodel velocity profiles
%**************************************************************************
%%                          MAIN ROUTINE                                  %
%**************************************************************************

% Initialize all global variables
Ui = mean(U_DNS(1,:,:),'all');       %use mean value of U-LES to approximate initial profile at all grid points on wall
Vi = mean(V_DNS(1,:,:),'all');
Pi = mean(p_DNS(1,:,:),'all');

initialize(nx,ny,dz,cfl,L_x,vonk,B,nu,iwm_tol,iwm_eps,Ui,Vi,Pi);

% Allocate size to wall-stress output arrays (only store values at one time
% instance)
txz(nx,ny,1)  = 0;
tyz(nx,ny,1)  = 0;

tic
% Loops over all time steps available in the input field arrays
for jt=1:jt_end
    
    fprintf('t=%5.4d\n',jt*dt)
    
    %LES solution (2 dimentional arrays)
    u(:,:)  = U_DNS(jt+1,:,:);
    v(:,:)  = V_DNS(jt+1,:,:);
    p(:,:)  = p_DNS(jt+1,:,:);
    w(:,:)  = W_DNS(jt+1,:,:);
    
    % time step seen by the wall model
    iwm_dt = dt;
    
    % Compute the wall stress
    
    % Update the convective term, turbulent diffusion term etc. on LHS of eqn
    % C22 in Yang et al 2015
    iwm_calc_lhs(nx,ny,dx,dy,u,v,p);
    
    % the subroutine to calculate wall stress
    iwm_calc_wallstress(nx,ny,iwm_tol,iwm_eps,MaxIter,nu,vonk,B,theta,jt*dt);
    
    % Imposing tauw_x, tauw_y in the LES solver as a boundary condition
    for iwm_i = 1:nx
        for iwm_j = 1:ny
            % wall stress, use the value calculated in iwm
            txz(iwm_i,iwm_j,jt) = iwm_tauwx(iwm_i,iwm_j);
            tyz(iwm_i,iwm_j,jt) = iwm_tauwy(iwm_i,iwm_j);
        end
    end
    
    % Collect some Stats
    tauw_tot(:,:,jt) = sqrt(txz(:,:,jt).^2 + tyz(:,:,jt).^2);
    tauw_wm(jt)      = mean(tauw_tot(:,:,jt),'all');
    eqmrate(jt)      = nnz(equil_flag)/numel(equil_flag);
    tR(jt)           = mean(iwm_tR(:,:),'all');
    
    % Flags for Diagnostics
    mat_eqm(jt,:)    = equil_flag(:,2);
    mat_div(jt,:)    = div_flag(:,2);
    mat_neg(jt,:)    = neg_flag(:,2);
    mat_maxiter(jt,:)= maxiter_flag(:,2);
    mat_A(jt,:)      = A_flag(:,2);
    
    divrate(jt,:)    = nnz(div_flag)/numel(div_flag);
    negrate(jt,:)    = nnz(neg_flag)/numel(neg_flag);
    Arate(jt,:)      = nnz(A_flag)/numel(A_flag);
    maxiterrate(jt,:)= nnz(maxiter_flag)/numel(maxiter_flag);
    
    % Mean velocity profile and instantaneous profile to get flow angle vs wall-normal
    % direction
    for iwm_i = 1:nx
        for iwm_j = 1:ny
            Vel = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)^2. +iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)^2.);
            utau = (iwm_utx(iwm_i,iwm_j)^4.+iwm_uty(iwm_i,iwm_j)^4. )^0.25;
            
            
            Dz = iwm_Dz(iwm_i,iwm_j);
            zz      = logspace(log10(11*nu/utau),log10(1),100);
            zz_shrt = logspace(log10(0.0001), log10(11*nu/utau), 50);
            
            if equil_flag(iwm_i,iwm_j)==0
                del_visc = nu*sqrt(iwm_utx(iwm_i,iwm_j)^2.+iwm_uty(iwm_i,iwm_j)^2.)/(utau^2.0);
                
                ulog(iwm_i,iwm_j,:)= utau*((iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel)*(log(zz/Dz)/vonk)  + iwm_Ax(iwm_i,iwm_j)*(zz/Dz) + iwm_Cx(iwm_i,iwm_j));
                uvisc(iwm_i,iwm_j,:)=iwm_utx(iwm_i,iwm_j)*zz_shrt/del_visc;
                vvisc(iwm_i,iwm_j,:)=iwm_uty(iwm_i,iwm_j)*zz_shrt/del_visc;
            else
                del_visc = nu/utau;
                
                ulog(iwm_i,iwm_j,:)= utau*(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel)*(log(zz/del_visc)/vonk + B);
                uvisc(iwm_i,iwm_j,:)=utau*(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel)*zz_shrt/del_visc;
                vvisc(iwm_i,iwm_j,:)=utau*(iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Vel)*zz_shrt/del_visc;
            end
            
        end
    end
    
    ulog_m(jt,:)  = mean(ulog ,[1 2]);
    uvisc_m(jt,:) = mean(uvisc,[1 2]);
    
end
toc

ulog_p = mean(ulog_m(10:end,:),1)/ut_av_DNS;
uvisc_p = mean(uvisc_m(10:end,:),1)/ut_av_DNS;

%% ************** PLOT RESULTS **************
%%%%%% Plot Profiles %%%%%%

semilogx(zz/(nu/ut_av_DNS),ulog_p,'r','LineWidth',2);
hold on
semilogx(zz_shrt/(nu/ut_av_DNS),uvisc_p,'b','LineWidth',2);
% semilogx(Dz/del_visc,u_log_Dz,'k+');      % at matching location
xlabel('y+');
ylabel('U+');
xlim([0.1 2000]);
ylim([0 26]);

%Import DNS data for plotting profile
ll     = xlsread('Channel_DNS_JHTDB_Retau1000/mean_data_jhu');   
zpll   = ll(:,2);
Upll   = ll(:,3);
semilogx(zpll,Upll,'k.');

legend('log-layer','viscous sublayer','mean-DNS (JHU Re_{\tau}=1000)','Location','NorthWest')
% legend('log-layer','viscous sublayer','At h^{+}_{wm}','mean-DNS (JHU Re_{\tau}=1000)','Location','NorthWest')


%%%%% Plot tauw ratio vs time %%%%%
mm       = xlsread('Channel_DNS_JHTDB_Retau1000/re-tau');
time_dns = mm(:,1);
ret_dns  = mm(:,2);
tw_DNS   = (nu*ret_dns).^2;

tt = 1:jt_end;
tratio_wm  = abs(tauw_wm)/ut_av_DNS^2;
tratio_DNS = tw_DNS/ut_av_DNS^2;

figure
plot(tt*dt,tratio_wm);
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
% subroutine Initialize
%**************************************************************************
% This subroutine initializes the wall model with plug flow
% conditions (assuming axial flow in x only)

function initialize(nx,ny,dz,cfl,L_x,vonk,B,nu,iwm_tol,iwm_eps,Ui,Vi,Pi)

global iwm_dirx iwm_diry iwm_DN iwm_Lu iwm_Luu iwm_Lv iwm_Lvv iwm_Luv ...
       iwm_utx iwm_uty iwm_tauwx iwm_tauwy iwm_flt_tagvel iwm_flt_tagvel_m iwm_flt_p  ...
       iwm_conv iwm_PrsGrad iwm_diff iwm_lhs iwm_flt_us ...
       iwm_tR iwm_Dz iwm_Ax iwm_Ay iwm_Cx iwm_Cy iwm_dt iwm_inte iwm_inte_m

% initial value for the x and y-velocity at first grid point
uinit = Ui;
vinit = Vi;

% at the height of the first grid point (h_wm)
Dz = dz/2.;

% cell height
iwm_Dz(1:nx,1:ny) = Dz; 

% filitered velocity at the first grid point in x, y directions
iwm_flt_tagvel  (1:nx,1:ny,iwm_dirx) = uinit;
iwm_flt_tagvel  (1:nx,1:ny,iwm_diry) = vinit;
iwm_flt_tagvel_m(1:nx,1:ny,iwm_dirx) = uinit;
iwm_flt_tagvel_m(1:nx,1:ny,iwm_diry) = vinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Note: We initialize utx and uty with standard log-law for equilibrium 

% set tolerance for
tol_eqm = iwm_tol;
eps_eqm = iwm_eps;

Vel = sqrt(uinit^2. + vinit^2.);
Uproj = uinit/Vel;
Vproj = vinit/Vel;

equilut = sign(Vel);
feqm = (equilut/vonk)*log(Dz*equilut/nu) + B*equilut - Vel;

iterm=0;
while abs(feqm) > tol_eqm && iterm <= 20
    ut_p = equilut + eps_eqm;
    feqm_ps = (ut_p/vonk)*log(Dz*ut_p/nu) + B*ut_p - Vel;
    a = (feqm_ps - feqm)/eps_eqm;
    equilut = equilut - feqm/a;
    
    feqm = (equilut/vonk)*log(Dz*equilut/nu) + B*equilut - Vel;
    iterm = iterm + 1;
end

equilutx = real(equilut*sqrt(Uproj));
equiluty = real(equilut*sqrt(Vproj));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

utau = (equilutx^4.+equiluty^4.)^0.25;
deli = 11*nu/utau;
delv = nu/utau;
rat  = deli/Dz;

% us in x, y directions
iwm_utx(1:nx,1:ny) = equilutx;
iwm_uty(1:nx,1:ny) = equiluty;

% wall stress in x, y directions
iwm_tauwx(1:nx,1:ny) = sign(equilutx)*equilutx^2.;
iwm_tauwy(1:nx,1:ny) = sign(equiluty)*equiluty^2.;

% pressure at first grid point
iwm_flt_p(1:nx,1:ny) = Pi;

% integrals of Lu, Lv, etc.
iwm_inte(1:nx,1:ny,iwm_Lu)  = utau*Uproj*(0.5*(deli^2.0)/delv + Dz*( B*(1-rat)...
                              - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(Dz/delv)) ));
        
iwm_inte(1:nx,1:ny,iwm_Lv)  = utau*Vproj*(0.5*(deli^2.0)/delv + Dz*( B*(1-rat)...
                              - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(Dz/delv)) ));
                          
iwm_inte(1:nx,1:ny,iwm_Luu) = (utau*Uproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
            Dz*((1.0/vonk)^2*((log(Dz/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
            +(2*B/vonk)*((log(Dz/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat))); 
        
iwm_inte(1:nx,1:ny,iwm_Lvv) = (utau*Vproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
            Dz*((1.0/vonk)^2*((log(Dz/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
            +(2*B/vonk)*((log(Dz/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));  
        
iwm_inte(1:nx,1:ny,iwm_Luv) = (utau^2)*Uproj*Vproj*(1.0/3.0*(deli^3.0/delv^2.0)+...
            Dz*((1.0/vonk^2)*((log(Dz/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
            +(2*B/vonk)*((log(Dz/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));  

iwm_inte_m(1:nx,1:ny,iwm_Lu)  = iwm_inte(1:nx,1:ny,iwm_Lu);
iwm_inte_m(1:nx,1:ny,iwm_Lv)  = iwm_inte(1:nx,1:ny,iwm_Lv);
iwm_inte_m(1:nx,1:ny,iwm_Luu) = iwm_inte(1:nx,1:ny,iwm_Luu);
iwm_inte_m(1:nx,1:ny,iwm_Lvv) = iwm_inte(1:nx,1:ny,iwm_Lvv);
iwm_inte_m(1:nx,1:ny,iwm_Luv) = iwm_inte(1:nx,1:ny,iwm_Luv);

% Each term in the integral equation and top/bottom derivatives
iwm_conv(1:nx,1:ny,1:iwm_DN)         = 0.;
iwm_PrsGrad(1:nx,1:ny,1:iwm_DN)      = 0.;
iwm_diff(1:nx,1:ny,1:iwm_DN)         = 0.;
iwm_lhs(1:nx,1:ny,iwm_dirx)          = -iwm_inte(1:nx,1:ny,iwm_Lu);     %ignore all terms multiplied by dt
iwm_lhs(1:nx,1:ny,iwm_diry)          = -iwm_inte(1:nx,1:ny,iwm_Lv);

% filtered friction velocity and the filtering time scale, tR<1
% Note: (cfl*L_x/nx/uinit) gives the LES time step
iwm_flt_us(1:nx,1:ny) = utau;
iwm_tR(1:nx,1:ny)     = (cfl*L_x/nx/uinit)/(Dz/vonk/equilutx); 

% linear correction to the log profile
iwm_Ax(1:nx,1:ny) = 0.;
iwm_Ay(1:nx,1:ny) = 0.;
iwm_Cx(1:nx,1:ny) = 0.;
iwm_Cy(1:nx,1:ny) = 0.;

% time step seen by the iwm
iwm_dt = (cfl*L_x/nx/uinit);

end

%**************************************************************************
%% subroutine iwm_calc_lhs()
%**************************************************************************
% This subroutine calculates the left hand side of the iwm system.

function iwm_calc_lhs(nx,ny,dx,dy,u,v,p)

global iwm_dirx iwm_diry iwm_Lu iwm_Luu iwm_Lv iwm_Lvv iwm_Luv ...
       iwm_flt_tagvel iwm_flt_tagvel_m iwm_flt_p  ...
       iwm_conv iwm_PrsGrad iwm_diff iwm_lhs ...
       iwm_tR iwm_Dz iwm_dt iwm_inte iwm_inte_m

% update u, v for the previous time step
for iwm_i = 1:nx
    for iwm_j = 1:ny
        iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx);
        iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry);
    end
end

% temporal filtering
for iwm_i=1:nx
    for iwm_j=1:ny
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx) =                               ...
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)*(1.-iwm_tR(iwm_i,iwm_j))    ...
            + u(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j);
        
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) =                               ...
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)*(1.-iwm_tR(iwm_i,iwm_j))    ...
            + v(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j);
        
        iwm_flt_p(iwm_i,iwm_j) = iwm_flt_p(iwm_i,iwm_j)                      ...
            * (1.-iwm_tR(iwm_i,iwm_j))                                       ...
            + p(iwm_i,iwm_j)*iwm_tR(iwm_i,iwm_j);
    end
end

% calculate LHS, calculation of the integrals is done from the last time step
% in the subroutine iwm_calc_wallstress, so is iwm_diff
for iwm_i = 1:nx
    for iwm_j = 1:ny
        
        % the convective term
        if iwm_i==1
            phip = iwm_inte(iwm_i+1,iwm_j,iwm_Luu);             
            phim = iwm_inte(nx,iwm_j,iwm_Luu);        %based on periodic condition with FVM
        elseif iwm_i==nx
            phip = iwm_inte(1,iwm_j,iwm_Luu);
            phim = iwm_inte(iwm_i-1,iwm_j,iwm_Luu);
        else
            phip = iwm_inte(iwm_i+1,iwm_j,iwm_Luu);
            phim = iwm_inte(iwm_i-1,iwm_j,iwm_Luu);
        end
        Luux = (phip-phim)/dx/2.;

        
        if iwm_j==1
            phip = iwm_inte(iwm_i,iwm_j+1,iwm_Luv);
            phim = iwm_inte(iwm_i,ny,iwm_Luv);
        elseif iwm_j==ny
            phip = iwm_inte(iwm_i,1,iwm_Luv);
            phim = iwm_inte(iwm_i,iwm_j-1,iwm_Luv);
        else
            phip = iwm_inte(iwm_i,iwm_j+1,iwm_Luv);
            phim = iwm_inte(iwm_i,iwm_j-1,iwm_Luv);
        end
        Luvy = (phip-phim)/dy/2.;

        
        if iwm_i==1
            phip = iwm_inte(iwm_i+1,iwm_j,iwm_Luv);
            phim = iwm_inte(nx,iwm_j,iwm_Luv);
        elseif iwm_i==nx
            phip = iwm_inte(1,iwm_j,iwm_Luv);
            phim = iwm_inte(iwm_i-1,iwm_j,iwm_Luv);
        else
            phip = iwm_inte(iwm_i+1,iwm_j,iwm_Luv);
            phim = iwm_inte(iwm_i-1,iwm_j,iwm_Luv);
        end
        Luvx = (phip-phim)/dx/2.;        
 
        
        if iwm_j==1
            phip = iwm_inte(iwm_i,iwm_j+1,iwm_Lvv);
            phim = iwm_inte(iwm_i,ny,iwm_Lvv);
        elseif iwm_j==ny
            phip = iwm_inte(iwm_i,1,iwm_Lvv);
            phim = iwm_inte(iwm_i,iwm_j-1,iwm_Lvv);
        else
            phip = iwm_inte(iwm_i,iwm_j+1,iwm_Lvv);
            phim = iwm_inte(iwm_i,iwm_j-1,iwm_Lvv);
        end
        Lvvy = (phip-phim)/dy/2.;

        
        if iwm_i==1
            phip = iwm_inte(iwm_i+1,iwm_j,iwm_Lu);
            phim = iwm_inte(nx,iwm_j,iwm_Lu);
        elseif iwm_i==nx
            phip = iwm_inte(1,iwm_j,iwm_Lu);
            phim = iwm_inte(iwm_i-1,iwm_j,iwm_Lu);
        else
            phip = iwm_inte(iwm_i+1,iwm_j,iwm_Lu);
            phim = iwm_inte(iwm_i-1,iwm_j,iwm_Lu);
        end
        Lux = (phip-phim)/dx/2.;
        
        
        if iwm_j==1
            phip = iwm_inte(iwm_i,iwm_j+1,iwm_Lv);
            phim = iwm_inte(iwm_i,ny,iwm_Lv);
        elseif iwm_j==ny
            phip = iwm_inte(iwm_i,1,iwm_Lv);
            phim = iwm_inte(iwm_i,iwm_j-1,iwm_Lv);
        else
            phip = iwm_inte(iwm_i,iwm_j+1,iwm_Lv);
            phim = iwm_inte(iwm_i,iwm_j-1,iwm_Lv);
        end
        Lvy = (phip-phim)/dy/2.;
        
        iwm_conv(iwm_i,iwm_j,iwm_dirx) = Luux + Luvy                               ...
            - iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_dirx)*(Lux+Lvy);
        iwm_conv(iwm_i,iwm_j,iwm_diry) = Luvx + Lvvy                               ...
            - iwm_flt_tagvel_m(iwm_i,iwm_j,iwm_diry)*(Lux+Lvy);
        
        % the pressure gradient term
        if iwm_i==1
            phip = iwm_flt_p(iwm_i+1,iwm_j);
            phim = iwm_flt_p(nx,iwm_j);
        elseif iwm_i==nx
            phip = iwm_flt_p(1,iwm_j);
            phim = iwm_flt_p(iwm_i-1,iwm_j);
        else
            phip = iwm_flt_p(iwm_i+1,iwm_j);
            phim = iwm_flt_p(iwm_i-1,iwm_j);
        end
        
        iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx) = (phip-phim)/dx/2.                ...
            * iwm_Dz(iwm_i,iwm_j);  %- 0.1*iwm_Dz(iwm_i,iwm_j);  %mean pressure gradient
        
        if iwm_j==1
            phip = iwm_flt_p(iwm_i,iwm_j+1);
            phim = iwm_flt_p(iwm_i,ny);
        elseif iwm_j==ny
            phip = iwm_flt_p(iwm_i,1);
            phim = iwm_flt_p(iwm_i,iwm_j-1);
        else
            phip = iwm_flt_p(iwm_i,iwm_j+1);
            phim = iwm_flt_p(iwm_i,iwm_j-1);
        end
        iwm_PrsGrad(iwm_i,iwm_j,iwm_diry) = (phip-phim)/dy/2.                ...
            * iwm_Dz(iwm_i,iwm_j);
        
        % the left hand side (this is f(y,t) in ODE: y'= f(y,t) , y = L )
        % this is the integrated momentum equation, except for the Lu term
        iwm_lhs(iwm_i,iwm_j,iwm_dirx) = -iwm_inte(iwm_i,iwm_j,iwm_Lu)              ...
            + iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_dirx)                              ...
            + iwm_PrsGrad(iwm_i,iwm_j,iwm_dirx)                                    ...
            - iwm_diff(iwm_i,iwm_j,iwm_dirx) );
        
        % this is the integrated momentum equation, except for the Lv term
        iwm_lhs(iwm_i,iwm_j,iwm_diry) = -iwm_inte(iwm_i,iwm_j,iwm_Lv)              ...
            + iwm_dt*( iwm_conv(iwm_i,iwm_j,iwm_diry)                              ...
            + iwm_PrsGrad(iwm_i,iwm_j,iwm_diry)                                    ...
            - iwm_diff(iwm_i,iwm_j,iwm_diry) );
    end
end

end

%**************************************************************************
%% subroutine iwm_calc_wallstress
%**************************************************************************

function iwm_calc_wallstress(nx,ny,iwm_tol,iwm_eps,MaxIter,nu,vonk,B,theta,t)

global iwm_dirx iwm_diry iwm_Lu iwm_Luu iwm_Lv iwm_Lvv iwm_Luv ...
       iwm_utx iwm_uty iwm_tauwx iwm_tauwy iwm_flt_tagvel ...
       iwm_diff iwm_lhs iwm_dudzT iwm_flt_us ...
       iwm_tR iwm_Dz iwm_Ax iwm_Ay iwm_Cx iwm_Cy iwm_dt iwm_inte iwm_inte_m ...
       equil_flag ...
       div_flag neg_flag A_flag maxiter_flag neg_flag_uty neg_flag_V deli_flag
       
   
for iwm_i=1:nx
for iwm_j=1:ny

% Note: We initialize utx and uty with standard log-law used in equilibrium flag 

    % set tolerance for 
    tol_eqm = iwm_tol;
    eps_eqm = iwm_eps;

    %total velocity
    Vel = sqrt(iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)^2.                  ...
        +iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)^2.);
    
    Uproj = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Vel;
    Vproj = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Vel;
        
    equilut = sign(Vel);
    
    feqm = (equilut/vonk)*log(iwm_Dz(iwm_i,iwm_j)*equilut/nu) + B*equilut - Vel;
    
    iterm=0;
    while abs(feqm) > tol_eqm && iterm <= 20
        ut_p = equilut + eps_eqm;
        feqm_ps = (ut_p/vonk)*log(iwm_Dz(iwm_i,iwm_j)*ut_p/nu) + B*ut_p - Vel;
        a = (feqm_ps - feqm)/eps_eqm;
        equilut = equilut - feqm/a;
        
        feqm = (equilut/vonk)*log(iwm_Dz(iwm_i,iwm_j)*equilut/nu) + B*equilut - Vel;
        iterm = iterm + 1;
    end   
    
    if iterm==21
        fprintf('ut_eqm did not converge in 20 iterations for i=%i j=%i\n at t=%f5.4',iwm_i,iwm_j,t);
    end

    deli_flag(iwm_i,iwm_j)=0;
    
    if abs(iwm_Dz(iwm_i,iwm_j)) <= abs(11*nu/equilut)
        equilut = sqrt(nu*Vel/iwm_Dz(iwm_i,iwm_j));
        deli_flag(iwm_i,iwm_j)=1;
    end
    
    equilutx = sign(Uproj)*equilut*sqrt(abs(Uproj));
    equiluty = sign(Vproj)*equilut*sqrt(abs(Vproj));
    
    %Seed values of utx & uty for interal WM calculation
    iwm_utx(iwm_i,iwm_j) = equilutx;
    iwm_uty(iwm_i,iwm_j) = equiluty;
    
    [fx,fy] = iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx), iwm_lhs(iwm_i,iwm_j,iwm_diry),   ...
        iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),    ...
        iwm_Dz(iwm_i,iwm_j), iwm_utx(iwm_i,iwm_j), iwm_uty(iwm_i,iwm_j), nu, vonk);
    
    iter = 0;
    equil_flag(iwm_i,iwm_j) = 0;
    div_flag(iwm_i,iwm_j) = 0;
    neg_flag(iwm_i,iwm_j) = 0;
    neg_flag_uty(iwm_i,iwm_j) = 0;
    neg_flag_V(iwm_i,iwm_j)   = 0;
    maxiter_flag(iwm_i,iwm_j) = 0;
    A_flag(iwm_i,iwm_j) = 0;

    
    if deli_flag(iwm_i,iwm_j)==0
    % use Newton method to solve the system
    while max(abs(fx),abs(fy)) > iwm_tol

        iwmutxP = iwm_utx(iwm_i,iwm_j)+iwm_eps;
        iwmutyP = iwm_uty(iwm_i,iwm_j);

        [fxp, fyp] = iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),iwm_lhs(iwm_i,iwm_j,iwm_diry), ...
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),...
            iwm_Dz(iwm_i,iwm_j), iwmutxP, iwmutyP, nu, vonk);

        a11 = (fxp-fx)/iwm_eps;
        a21 = (fyp-fy)/iwm_eps;

        iwmutxP = iwm_utx(iwm_i,iwm_j);
        iwmutyP = iwm_uty(iwm_i,iwm_j)+iwm_eps;

        [fxp, fyp] = iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx),iwm_lhs(iwm_i,iwm_j,iwm_diry), ...
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),    ...
            iwm_Dz(iwm_i,iwm_j), iwmutxP, iwmutyP, nu, vonk);

        a12 = (fxp-fx)/iwm_eps;
        a22 = (fyp-fy)/iwm_eps;
        
        iwm_utx(iwm_i,iwm_j) = iwm_utx(iwm_i,iwm_j) - 0.50*( a22*fx-a12*fy)/(a11*a22-a12*a21);
        iwm_uty(iwm_i,iwm_j) = iwm_uty(iwm_i,iwm_j) - 0.50*(-a21*fx+a11*fy)/(a11*a22-a12*a21);
        
        % infinity check
        if (abs(iwm_utx(iwm_i,iwm_j)) > 10^4.0 || abs(iwm_uty(iwm_i,iwm_j))> 10^4.0)
            fprintf('divergence in utx/uty at i=%i j=%i\n',iwm_i,iwm_j);            
            equil_flag(iwm_i,iwm_j) = 1;
            div_flag(iwm_i,iwm_j)   = 1;
            break
        end
        
        if iwm_utx(iwm_i,iwm_j) < 0
%             fprintf('utx has a negative value at i=%i j=%i\n',iwm_i,iwm_j);
            neg_flag(iwm_i,iwm_j)   = 1;
        end
        
        if iwm_uty(iwm_i,iwm_j) < 0
            neg_flag_uty(iwm_i,iwm_j)   = 1;
        end 
        
        if iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry) < 0
            neg_flag_V(iwm_i,iwm_j)   = 1;
        end
        
        [fx, fy] = iwm_slv(iwm_lhs(iwm_i,iwm_j,iwm_dirx), iwm_lhs(iwm_i,iwm_j,iwm_diry),   ...
            iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx), iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry),    ...
            iwm_Dz(iwm_i,iwm_j), iwm_utx(iwm_i,iwm_j), iwm_uty(iwm_i,iwm_j), nu, vonk);
 
        iter = iter+1;

        % maximum iteration reached
        if (iter>MaxIter)
            equil_flag(iwm_i,iwm_j) = 1;
            maxiter_flag(iwm_i,iwm_j) = 1;
            break
        end
    end
    end
    
    %SWITCH TO EQUILIBRIUM MODEL 
    if (equil_flag(iwm_i,iwm_j)==1 || deli_flag(iwm_i,iwm_j)==1)
        iwm_utx(iwm_i,iwm_j) = equilutx;
        iwm_uty(iwm_i,iwm_j) = equiluty;
    end
    
    %%%% At this point utx and uty have been calculated. After this step, we simply
    %calculate other parameters using their relations with utx, uty and utau
    
    % calculate the friciton velocity
    utau = (iwm_utx(iwm_i,iwm_j)^4.+iwm_uty(iwm_i,iwm_j)^4. )^0.25;
        
    % eq. C25 in Yang et al. 2015
    deli = min( 11*nu/utau , iwm_Dz(iwm_i,iwm_j) );    
    rat = (deli/iwm_Dz(iwm_i,iwm_j));
    
    %calculate Ax, Ay, Cx, Cy
    if(equil_flag(iwm_i,iwm_j)==1  || deli_flag(iwm_i,iwm_j)==1)
        %By definition of single profile characterized by utau for eqm flag
        delv = nu/utau;
        Ax = 0.;
        Ay = 0.;
        Cx = 0.;
        Cy = 0.;
    else
        % eq. C26 in Yang et al. 2015 with Hayat & Park 2021 modification

        delv = nu/utau;
        
        if deli == 11*nu/utau
            % eq. C27 in Yang et al. 2015 with Hayat & Park 2021 modification
            Ax = (iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/utau + Uproj/vonk*log(rat) - (deli/delv)*sign(iwm_utx(iwm_i,iwm_j))*(iwm_utx(iwm_i,iwm_j)/utau)^2.0 ) / ((1.0-rat));
            Ay = (iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/utau + Vproj/vonk*log(rat) - (deli/delv)*sign(iwm_uty(iwm_i,iwm_j))*(iwm_uty(iwm_i,iwm_j)/utau)^2.0 ) / ((1.0-rat));
            Cx = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/utau - Ax;
            Cy = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/utau - Ay;
            
        elseif deli == iwm_Dz(iwm_i,iwm_j)
            Ax = 0.;
            Ay = 0.;
            Cx = 0.;
            Cy = 0.;
            
        end
    end

    % check for excessive linear term correction
    % Revert to EQWM if values exceed 2.0
    if (abs(Ax)>2|| abs(Ay)>2)
        equil_flag(iwm_i,iwm_j) = 1;
        A_flag(iwm_i,iwm_j) = 1;
        iwm_utx(iwm_i,iwm_j) = equilutx;
        iwm_uty(iwm_i,iwm_j) = equiluty;
        Ax = 0.;
        Ay = 0.;
        Cx = 0.;
        Cy = 0.;
        
        utau = (iwm_utx(iwm_i,iwm_j)^4.+iwm_uty(iwm_i,iwm_j)^4. )^0.25;
        deli = 11*nu/utau;
        
        if deli_flag(iwm_i,iwm_j)==1
        deli = iwm_Dz(iwm_i,iwm_j);    
        end
        
        delv = nu/utau;
        rat = (deli/iwm_Dz(iwm_i,iwm_j));        
    end
    
    % store the linear correction
    iwm_Ax(iwm_i,iwm_j) = Ax;
    iwm_Ay(iwm_i,iwm_j) = Ay;
    iwm_Cx(iwm_i,iwm_j) = Cx;
    iwm_Cy(iwm_i,iwm_j) = Cy;
    
    % update integral for last time step
    iwm_inte_m(iwm_i,iwm_j,iwm_Lu ) = iwm_inte(iwm_i,iwm_j,iwm_Lu );
    iwm_inte_m(iwm_i,iwm_j,iwm_Lv ) = iwm_inte(iwm_i,iwm_j,iwm_Lv );
    iwm_inte_m(iwm_i,iwm_j,iwm_Luv) = iwm_inte(iwm_i,iwm_j,iwm_Luv);
    iwm_inte_m(iwm_i,iwm_j,iwm_Luu) = iwm_inte(iwm_i,iwm_j,iwm_Luu);
    iwm_inte_m(iwm_i,iwm_j,iwm_Lvv) = iwm_inte(iwm_i,iwm_j,iwm_Lvv);

    % calculate the required integrals

    Dz = iwm_Dz(iwm_i,iwm_j);     %To make following equations more legible
    utau_x = iwm_utx(iwm_i,iwm_j);
    utau_y = iwm_uty(iwm_i,iwm_j);
    
    if (equil_flag(iwm_i,iwm_j)==1 || deli_flag(iwm_i,iwm_j)==1)
        if deli_flag(iwm_i,iwm_j)==0
            % Lu
            iwm_inte(iwm_i,iwm_j,iwm_Lu) = utau*Uproj*(0.5*(deli^2.0)/delv + ...
                Dz*( B*(1-rat) - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(Dz/delv)) ));
            % Lv
            iwm_inte(iwm_i,iwm_j,iwm_Lv) = utau*Vproj*(0.5*(deli^2.0)/delv + ...
                Dz*( B*(1-rat) - 1.0/vonk*(1-rat+rat*log(deli/delv)-log(Dz/delv)) ));
            % Luu
            iwm_inte(iwm_i,iwm_j,iwm_Luu) = (utau*Uproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
                Dz*((1.0/vonk)^2*((log(Dz/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
                +(2*B/vonk)*((log(Dz/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));
            % Lvv
            iwm_inte(iwm_i,iwm_j,iwm_Lvv) = (utau*Vproj)^2 * (1.0/3.0*(deli^3.0/delv^2.0)+...
                Dz*((1.0/vonk)^2*((log(Dz/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
                +(2*B/vonk)*((log(Dz/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));
            % Luv
            iwm_inte(iwm_i,iwm_j,iwm_Luv) = (utau^2)*Uproj*Vproj*(1.0/3.0*(deli^3.0/delv^2.0)+...
                Dz*((1.0/vonk^2)*((log(Dz/delv) - 1)^2.0 -rat*(log(deli/delv) - 1)^2.0 + 1-rat) ...
                +(2*B/vonk)*((log(Dz/delv) - 1) -rat*(log(deli/delv) - 1)) + B^2.0*(1-rat)));
        else
            % Lu:  Eq. C19 in Yang et al. 2015 with Hayat & Park 2021 modification
            iwm_inte(iwm_i,iwm_j,iwm_Lu) = 0.5*sign(utau_x)*((utau_x^2.0)/utau)*(deli^2.0)/delv;
            % Lv
            iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0.5*sign(utau_y)*((utau_y^2.0)/utau)*(deli^2.0)/delv;
            
            % Luu: Eq. C20 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1.0/3.0*(utau_x^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            % Lvv
            iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1.0/3.0*(utau_y^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            
            % Luv: Eq. C21 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1.0/3.0*sign(utau_x)*sign(utau_y)*(utau_x*utau_y/utau)^2.0*(deli^3.0/delv^2.0);
        end
        
    else
        
        if deli == 11*nu/utau
            
            % Lu:  Eq. C19 in Yang et al. 2015 with Hayat & Park 2021 modification
            iwm_inte(iwm_i,iwm_j,iwm_Lu) = 0.5*sign(utau_x)*((utau_x^2.0)/utau)*(deli^2.0)/delv + ...
                utau*Dz*(0.5*Ax*(1-rat^2.0) + Cx*(1-rat) - Uproj/vonk*(1-rat+rat*log(rat)));
            % Lv
            iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0.5*sign(utau_y)*((utau_y^2.0)/utau)*(deli^2.0)/delv + ...
                utau*Dz*(0.5*Ay*(1-rat^2.0) + Cy*(1-rat) - Vproj/vonk*(1-rat+rat*log(rat)));
            
            % Luu: Eq. C20 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1.0/3.0*(utau_x^4.0/utau^2.0)*(deli^3.0/delv^2.0)+...
                utau^2.0*Dz*(-Ax*(Uproj/vonk)*rat^2.0*log(rat) + Ax*(Cx-Uproj/vonk/2.0) ...
                *(1-rat^2.0) + Ax^2.0/3.0*(1-rat^3.0) + (Cx - Uproj/vonk)^2.0   ...
                -rat*(Cx-Uproj/vonk+Uproj/vonk*log(rat))^2.0 + (Uproj/vonk)^2.0*(1-rat));
            % Lvv
            iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1.0/3.0*(utau_y^4.0/utau^2.0)*(deli^3.0/delv^2.0)+...
                utau^2.0*Dz*(-Ay*(Vproj/vonk)*rat^2.0*log(rat) + Ay*(Cy-Vproj/vonk/2.0) ...
                *(1-rat^2.0) + Ay^2.0/3.0*(1-rat^3.0) + (Cy - Vproj/vonk)^2.0   ...
                -rat*(Cy-Vproj/vonk+Vproj/vonk*log(rat))^2.0 + (Vproj/vonk)^2.0*(1-rat));
            
            % Luv: Eq. C21 in Yang et al 2015 with Hayat & Park 2021 modification
            iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1.0/3.0*sign(utau_x)*sign(utau_y)*(utau_x*utau_y/utau)^2.0*(deli^3.0/delv^2.0)+...
                utau^2.0*Dz*(-1.0/vonk*(Ax*Vproj+Ay*Uproj)*(0.25-0.25*rat^2.0+0.5*rat^2.0*log(rat))...       %CHECK last rat^2
                -1/vonk*(Cx*Vproj+Cy*Uproj)*(1-rat+rat*log(rat))-Uproj*Vproj/(vonk^2.0)*(rat-2+rat*(log(rat)-1)^2.0)...
                +(1.0/3.0)*Ax*Ay*(1-rat^3.0)+0.5*(Ax*Cy+Ay*Cx)*(1-rat^2.0)+Cx*Cy*(1-rat));
            
        elseif deli == Dz
            
            % Lu:  
            iwm_inte(iwm_i,iwm_j,iwm_Lu) = 0.5*sign(utau_x)*((utau_x^2.0)/utau)*(deli^2.0)/delv;
            % Lv
            iwm_inte(iwm_i,iwm_j,iwm_Lv) = 0.5*sign(utau_y)*((utau_y^2.0)/utau)*(deli^2.0)/delv;
            
            % Luu:
            iwm_inte(iwm_i,iwm_j,iwm_Luu) = 1.0/3.0*(utau_x^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            % Lvv
            iwm_inte(iwm_i,iwm_j,iwm_Lvv) = 1.0/3.0*(utau_y^4.0/utau^2.0)*(deli^3.0/delv^2.0);
            
            % Luv: 
            iwm_inte(iwm_i,iwm_j,iwm_Luv) = 1.0/3.0*sign(utau_x)*sign(utau_y)*(utau_x*utau_y/utau)^2.0*(deli^3.0/delv^2.0);
        end
    end

    if  deli_flag(iwm_i,iwm_j)==1
        
        % calculate top derivatives
        iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Dz;
        iwm_dudzT(iwm_i,iwm_j,iwm_diry) = iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Dz;
        
        % calculte the turbulent diffusion term
        dVelzT = Uproj*iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx)/Dz + Vproj*iwm_flt_tagvel(iwm_i,iwm_j,iwm_diry)/Dz;        

    else
        
        % calculate top derivatives
        % From Eq. C12
        iwm_dudzT(iwm_i,iwm_j,iwm_dirx) = (Uproj/vonk+Ax)*utau/Dz;
        iwm_dudzT(iwm_i,iwm_j,iwm_diry) = (Vproj/vonk+Ay)*utau/Dz;
        
        % calculte the turbulent diffusion term
        % Eq. C14
        %Note: this reduces to utau/(vonk*Dz) for eqm flag
        dVelzT = (Uproj*(Uproj/vonk+Ax) + Vproj*(Vproj/vonk+Ay))*utau/Dz;        
    end
               
    % Eq. C13, the eddy viscosity 'nut' at y=hwm is:
    nut_hwm = (vonk*Dz)^2.0*dVelzT;
    
    % calculate the wall stress
    %Note: Unlike prescribed roughness case, here we have the same form for
    %tauwx for both equilibrium and non-equilibrium
    
    iwm_tauwx(iwm_i,iwm_j) = sign(iwm_utx(iwm_i,iwm_j))*iwm_utx(iwm_i,iwm_j)^2.0;
    iwm_tauwy(iwm_i,iwm_j) = sign(iwm_uty(iwm_i,iwm_j))*iwm_uty(iwm_i,iwm_j)^2.0;
    
    % shear stress difference between top and bottom
    iwm_diff(iwm_i,iwm_j,iwm_dirx) = (nut_hwm+nu)*iwm_dudzT(iwm_i,iwm_j,...
        iwm_dirx) - iwm_tauwx(iwm_i,iwm_j);
    
    iwm_diff(iwm_i,iwm_j,iwm_diry) = (nut_hwm+nu)*iwm_dudzT(iwm_i,iwm_j,...
        iwm_diry) - iwm_tauwy(iwm_i,iwm_j);
    
    % TIMESCALE: the filtered friction velocity used for filtering time scale
    iwm_flt_us(iwm_i,iwm_i) = iwm_tR(iwm_i,iwm_j)*utau + (1.-iwm_tR(iwm_i,iwm_j))*iwm_flt_us(iwm_i,iwm_j);
    
    % update the filtering time scale E
    % Eq. 26
    
    iwm_tR(iwm_i,iwm_j) = iwm_dt/(theta*Dz/iwm_flt_us(iwm_i,iwm_j)/vonk);
    
    % filtering time scale can only be larger than the time step,
    % if not, then just use the instantaneous flow field to do the model
    if (iwm_tR(iwm_i,iwm_j) > 1.)
        iwm_tR(iwm_i,iwm_j) = 1.;
    end
end
end

end
%*******************************************************************************
%% subroutine iwm_slv(lhsx,lhsy,iwm_flt_tagvel(iwm_i,iwm_j,iwm_dirx),Uy,Dz,deli,utx,uty,fx,fy)
%*******************************************************************************

function [fx,fy]= iwm_slv(lhsx,lhsy,Ux,Uy,Dz,utx,uty,nu,vonk)

%U_tot
Vel = sqrt(Ux^2.0+Uy^2.0);

%Eqn C30
utau = (utx^4.0 + uty^4.0)^(0.25);

%Eqn C26 (viscous scale)
delv = nu/utau;

%Eqn C25 
%(Interface of log and viscous layer chosen as min of y+=11 or h_wm)
deli = min(11*nu/utau,Dz);
% deli = 11*nu/utau;             % 11*delv; 

%Find parameters of log profile only if deli=11*nu/utau, otherwise linear
%profile in whole inner layer
if deli == 11*nu/utau
    %Eqn C27
    Ax = (Ux/utau + (Ux/Vel)/vonk*log(deli/Dz) - (deli/delv)*sign(utx)*(utx/utau)^2.0) / ((1.0-deli/Dz));
    Ay = (Uy/utau + (Uy/Vel)/vonk*log(deli/Dz) - (deli/delv)*sign(uty)*(uty/utau)^2.0) / ((1.0-deli/Dz));
    
    %Eqn C28
    Cx = Ux/utau - Ax;
    Cy = Uy/utau - Ay;
    
    %LHS of Eqn C23 and C24
    %define ratio of deli to h_wm
    rat = (deli/Dz);
    inteLu = 0.5*sign(utx)*(utx^2.0/utau)*(deli^2.0)/delv + utau*Dz*( 0.5*Ax*(1-rat^2.0) + Cx*(1-rat) - (Ux/Vel)/vonk*(1-rat+rat*log(rat)) );
    inteLv = 0.5*sign(uty)*(uty^2.0/utau)*(deli^2.0)/delv + utau*Dz*( 0.5*Ay*(1-rat^2.0) + Cy*(1-rat) - (Uy/Vel)/vonk*(1-rat+rat*log(rat)) );
    
elseif deli == Dz
    
    inteLu = 0.5*sign(utx)*(utx^2.0/utau)*(deli^2.0)/delv;
    inteLv = 0.5*sign(uty)*(uty^2.0/utau)*(deli^2.0)/delv; 
    
end

fx = inteLu+lhsx;
fy = inteLv+lhsy; 

end
