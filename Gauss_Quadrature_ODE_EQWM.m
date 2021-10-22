%********  Gauss-Quadrature based ODE Equilibrium Wall Model ************%
%                                                                        %
% Reference:                                                             %
% [1] "Numerical implementation of efficient grid-free integral          %
%     wallmodels in unstructured-grid LES solvers" (Hayat & Park, 2021)  %
%                                                                        %
% Written by: Imran Hayat - 10/21/2021                                   %
%                                                                        %
%                                                                        %
%************************************************************************%

clear all
clc

%% Main inputs

% Number of quadrature points                     
% For GLL, n>=(N+3)/2 for exact integral of polynomial order N
n=[30]; 

% Choose Re_tau at which DNS mean profile data will be imported 
% (source: https://turbulence.oden.utexas.edu/)
Re = [550 950 2000 4200 5200];  

%% Secondary inputs
iter = 10;           % For Secant method
tol  = 1e-8;         % For Secant method
                   
%% Parameters
K     = 0.41;
Aplus = 26;
B     = 5.3;
rho   = 1.0;

%Boundary points and Boundary conditions
y_0   = 0;              % at wall
y_hwm = 0.1;            % at matching height    
u0    = 0;              % no-slip
              
% Note: u_hwm is defined below inside the loop for Retau

for r = 1:length(Re)

if Re(r)==180
    index=1;
elseif Re(r)==550
    index=2;
elseif Re(r)==950
    index=3;
elseif Re(r)==2000
    index=4;
elseif Re(r)==4200
    index=5;    
elseif Re(r)==5200
    index=6;    
end

% read in DNS y+ vs u+ data to obtain u_hwm at y_hwm
filename = sprintf('./DNS_data/Austin_Retau%i.dat',Re(r));
data  = load(filename);
y     = data(:, 1);
yplus = data(:, 2);
uplus = data(:, 3);

% Read in corresponding flow variables to dimensionalize the data
% ("Austin_DNS_param.xlsx" contains kinematic viscosity and utau
% information at different Re_tau obtained from UT-Austin DNS database.)

param = xlsread('./DNS_data/Austin_DNS_param.xlsx');
utau_nominal = param(index,3);         
Retau        = param(index,1);

mu_lam =  utau_nominal/Retau;            
nu     =  mu_lam/rho;
u      =  uplus * utau_nominal;
y      =  yplus * nu / utau_nominal;

% interpolate u at chosen y_hwm
u_hwm = interp1(y, u, y_hwm);
clear yplus uplus data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ********* CODE FOR GAUSS-QUADRATURE ODE WALL MODEL STARTS HERE *******%%

%% Guess values for utau

% Get the first guess value from linear profile
ut1 = sqrt(u_hwm*nu/y_hwm);         

% Get the second guess value from the log profile
% For this guess value, solve iteratively using Newton-Raphson
tol_NR  = 1e-5;
eps_eqm = 1e-5;
ut2 = ut1;

% Define log profile function
f_ut2 = @(ut) ( (ut/K)*log(y_hwm*ut/nu) + B*ut - u_hwm   );

ut2 = newton_raphson(ut2,f_ut2,tol_NR,eps_eqm);

tic
%% Gauss-Labatto-Legendre (GLL or LGL)
% Returns quadrature pts x and weights w in the interval [-1,1] 

% (Note: The fuction lglnodes.m was obtained from the following repository:
%
% https://github.com/ClimateGlobalChange/SpectralElementMethodAnalysis...
% /blob/master/lglnodes.m
%
% Alternately, one can use standard Gauss-Labatto-Legendre tables to find 
% coordinates and corresponding weights for a specified n.)

%  [x,w,TT2] = lglnodes(y_0,y_hwm,n);     
 [x_lgl,w_lgl,TT2] = lglnodes(n); 


% Transform to original coordinates [y_0,y_hwm] from legendre interval [-1,1]
% We use exponential transformation which clusters more points close to the
% wall.

% linear transform (no clustering, eqn (3.2) in ref [1]) :
% x = (y_hwm+y_0)/2 + (y_hwm-y_0).*x_lgl/2;
% w = (y_hwm-y_0)/2 * w_lgl;

% nonlinear transform (allows clustering, eqn (3.4) in ref [1]):
x = y_hwm*(exp(x_lgl+1)-1)/(exp(2)-1);
w = y_hwm*w_lgl.*exp(x_lgl+1)/(exp(2)-1);

%% Integral Evaluation using Gauss Quadrature

% We use the expression for Integrand defined in eqn (2.14) in ref [1].

% Define Integrand
lm_p = @(y,utau) ( (K*utau/nu)*y*(1-exp(-utau*y/(nu*Aplus)))  );
f    = @(y,utau) ( (2*utau^2/nu)/(1+sqrt(1+4*lm_p(y,utau)^2)) ); %Integrand

%% Shooting Method (Secant)
% Use secant method to converge the value of integral

[utnew,iter_gq,time_elapsed] = quad_secant(x,w,u_hwm,f,ut1,ut2,tol,iter);

utau(r)  = utnew;
err(r)   = 100*abs(utau_nominal^2-utnew^2)/utau_nominal^2;


%% Plot velocity profiles

xplot = linspace(y_0,y_hwm,1000);

for k=1:length(xplot)
    [zint,wwint] = lglnodes(n);
    
    xint = xplot(k)*(exp(zint+1)-1)/(exp(2)-1);
    wint = xplot(k)*wwint.*exp(zint+1)/(exp(2)-1);  
    
    u_wm(k) = 0;
    for j=1:length(xint)
        u_wm(k) = u_wm(k) + f(xint(j),utnew) * wint(j);
    end
    
end


figure(r+100)    
semilogx(y/(nu/utau_nominal),u/utau_nominal, 'k.'); 
hold on;
semilogx(xplot/(nu/utau_nominal), u_wm/utau_nominal,'linewidth', 1.5); 
xlabel('y^{+}')
ylabel('U^{+}')

legendCell = strcat('n=',string(num2cell(n)));
addDNS=['DNS' legendCell];
legend(addDNS,'Location','NorthWest')
ttl=sprintf('Re_{t}=%i',Re(r));
title(ttl)
xlim([0.1 5000])

end

%% Plot error in tau_w
figure(200)
hold on
plot(Re, err','-o', 'linewidth', 1.5);
xlabel('Re_{\tau}')
ylabel('% Error in \tau_{w}')
legend(legendCell)
ylim([0 5])


%************************ END OF GQWM CODE *******************************%




%*************************************************************************%
%                           Subroutines                                   %
%*************************************************************************%

% Newton-Raphson Method:
function ut = newton_raphson(u,f,tol_NR,eps_eqm)

func = f(u);
iterm = 0;

while abs(func) > tol_NR && iterm <= 5
    up     = u + eps_eqm;
    func_p = f(up);
    slope  = (func_p - func)/eps_eqm;
    u      = u - func/slope;    
    func   = f(u);
    iterm  = iterm + 1;
end

ut = u;

end


% Secant method for quadrature formula 
function [utnew,iter_completed,T] = quad_secant(x,w,u_hwm,f,ut1,ut2,tol,iter)

tic
% Initial Solves
I1=0;
I2=0;

for j=1:length(x)
    I1 = I1 + f(x(j),ut1) * w(j);
    I2 = I2 + f(x(j),ut2) * w(j);
end

for i = 1:iter
    Inew=0;
    utnew = ut1 + (ut1-ut2)/(I1-I2)*(u_hwm-I1);
        
    for j=1:length(x)
        Inew = Inew + f(x(j),utnew) * w(j);
    end

    if ( abs(Inew-u_hwm) <= tol )
        break;
    end
    
    ut2 = ut1;
    ut1 = utnew;
    I2 = I1;
    I1 = Inew;
end
T=toc;
iter_completed = i;


if i==iter
    disp('max iteration reached in Secant method')
end 

end





% Funtion to generate GLL zeroes and weights
function [x,w,T]=lglnodes(N)
N=N-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods. 
%
% Reference on LGL nodes and weights: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
% Modified by Paul Ullrich - 10/12/2017
% Changed input to 
%
% Copyright (c) 2009, Greg von Winckel 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% * Redistributions of source code must retain the above copyright 
% notice, this list of conditions and the following disclaimer. 
% * Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in 
% the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
             
end

w=2./(N*N1*P(:,N1).^2);

x=flip(x);
w=flip(w);
T=toc;

end