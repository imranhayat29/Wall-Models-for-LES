%************************************************************************%
%                                                                        %
%       Traditional Finite-Volume based ODE Equilibrium Wall Model       %
%                                                                        %
%************************************************************************%

close all;
clc; clear all;
warning('off','all');

%% Wall model set up
%Flow conditions:
Re_tau = 2000;
rho    = 1.0;
mu_lam = 1.0/Re_tau;
nu_lam = mu_lam/rho;

%WM grid info:
ncv_wm         = 30;
stretch_factor = 1.1;

%Matching location:
h_wm_arr = 0.1;
h_wm_arr = h_wm_arr';

% read in DNS y+ vs u+ profile
filename = sprintf('./Channel_DNS_data/Austin_Retau%i.dat',Re_tau);  
data = load(filename);

% placeholder to store wall shear stress from the wall model. 
tauwall_topbc_arr = h_wm_arr*0;

% DNS data below 
y     = data(:, 1);
yplus = data(:, 2);
uplus = data(:, 3);
utau_nominal = 1.0;
u = uplus * utau_nominal;
y = yplus * nu_lam / utau_nominal;
     
%%
% Dirichlet top bc at y=h_wm.

for k=1:length(h_wm_arr)
	h_wm = h_wm_arr(k);
	u_les = interp1(y, u, h_wm );
    
	[u_wm_topbc, y_cv_topbc,  tauwall_topbc, iter] = solveUOnly_eqwm_Dirichlet_topBC(u_les, mu_lam, rho, h_wm, ncv_wm, stretch_factor); 
    
    tauwall_topbc_arr(k) = tauwall_topbc;
    utau(k)  = sqrt(tauwall_topbc_arr(k));
end

%% Plotting velocity profile
% figure();
% plot(y, u, 'k.'); hold on;
% plot(y_cv_topbc, u_wm_topbc, 'linewidth', 2); hold on;

figure();
semilogx(y/(nu_lam / utau_nominal), u/utau_nominal, 'k.'); hold on;
semilogx(y_cv_topbc/(nu_lam / utau_nominal), u_wm_topbc/utau_nominal, 'linewidth', 2); hold on;
xlim([0.01 5000])
xlabel('y^{+}')
ylabel('U^{+}')

error = 100*abs(utau.^2 - utau_nominal.^2)./(utau_nominal.^2);
fprintf('Error in Tau_wall=%4.2f %% \n',error);


%%
%*************************************************************************%
%                           Functions                                     %
%*************************************************************************%

% function to solve ODE Equilibrium wall model

function [u_wm, y_cv, tau_wall, iter] = solveUOnly_eqwm_Dirichlet_topBC( umag_les, mu_lam, rho, h_wm, ncv, stretch_factor )

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% h_wm is the y-coord. of the highest CV ( h_wm == max(y_cv) ), 
	% where the reference state (LES or DNS) is to be enforced. 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	kappa = 0.41;
	Aplus = 17.0;
    
	tol = 1e-8;
	max_iter = 2000;

	% reference laminar solution
	tau_wall_lam = mu_lam*umag_les/h_wm;

	% build implicit WM grid for ODE solve 
	[y_cv, y_fa] = buildGrid_eqwm_dirichlet_topbc( ncv, stretch_factor, h_wm);
	nfa = size(y_fa)*[1;0];

	iter = 0;
	done = 0;
	tau_wall = 0.7;

	while (done == 0)

		tau_wall_prev = tau_wall;

		% populate total viscosity at face 
		mu_fa = mu_lam * ones(nfa, 1);
		for ifa=2:nfa
			nu = mu_lam / rho;
			utau = sqrt(tau_wall / rho);
			D = 1.0-exp(-1.0*y_fa(ifa)*utau/nu/Aplus);
			D = D*D;
			mut = rho * kappa * y_fa(ifa) * utau * D;
			mu_fa(ifa) = mu_fa(ifa) + mut;
		end
		assert( mu_fa(1) == mu_lam );
		
		% construct momentum matrix system
		% solution vector is u_wm at y_cv. 
		A = zeros(ncv,ncv);
		b = zeros(ncv, 1);
		% top bc .. 
		A(ncv, ncv) = 1.0;		b(ncv, 1) = umag_les;
		% internal cvs 
		for icv = 2:ncv-1
			superdiag = mu_fa(icv+1)/( y_cv(icv+1) - y_cv(icv  ) );
			subdiag   = mu_fa(icv  )/( y_cv(icv  ) - y_cv(icv-1) );
			diag = -1.0* (superdiag+subdiag);
			A(icv,icv) = diag;
			A(icv, icv+1) = superdiag;
			A(icv, icv-1) = subdiag;
			b(icv) = 0.0;
		end
		% wall bc (enforced in eqn for icv=1, wall CV). 
		b(1) = 0.0;
		superdiag = mu_fa(2)/( y_cv(2)-y_cv(1) );
		diag = -1.0* ( superdiag + mu_fa(1)/( y_cv(1)-y_fa(1) ) );
		A(1,1) = diag;
		A(1,2) = superdiag; 

		% invert the momentum system 
		u_wm = A\b;
		% update tau_wall
		tau_wall = mu_fa(1) * ( u_wm(1) - 0.0)/(y_cv(1)-y_fa(1));

		if ( abs( (tau_wall - tau_wall_prev)/tau_wall_lam) < tol   )
			y1_plus = y_cv(1)/nu*sqrt(tau_wall/rho);
			assert(y1_plus < 1.0);
			done = 1;
		end
		
		if (done == 0 )
			iter = iter + 1;
		end

	end

end



% function to build grid for wall model

function [ y_cv, y_fa ] = buildGrid_eqwm_dirichlet_topbc(ncv, stretch_factor, h_wm)

    % velocity is stored at the cell center (y). 
		% In the Dirichlet top bc formulation, h_wm is the centroidal y-coord. of the last cell, 
		% where LES condition is to be enforced. 

    assert( h_wm > 0);
    assert( ncv > 0);
    assert( stretch_factor >= 1.0); 
    
		nfa = ncv + 1;
    y_cv = zeros(ncv, 1);
    y_fa = zeros(nfa, 1);
    
    y_fa(1) = 0;
    tmp_inc = 1.0;
    
    for ifa = 2:nfa
			y_fa(ifa) = y_fa(ifa-1) + tmp_inc;
			y_cv(ifa-1) = 0.5* ( y_fa(ifa) + y_fa(ifa-1) ); 
			tmp_inc = tmp_inc*stretch_factor; 
		end

		% rescale grid so that y_cv(max) = h_wm
		% LES condition is enforced at the centroid of the last cell. 
		y_max = max(y_cv);
		for ifa = 1:nfa
			y_fa(ifa) = y_fa(ifa)/y_max*h_wm;
		end
		for icv = 1:ncv
			y_cv(icv) = y_cv(icv)/y_max*h_wm;
		end
    

end

