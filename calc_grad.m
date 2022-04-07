function [phi_x,phi_z] = calc_grad(nx,nz,dx,dz,phi)
% Computes the gradients in streamwise (x) and spanwise (z) directions
% using finite differencing on a cartesian grid with periodic boundary 
% conditions assumed.

% nx    = number of grid points in x
% dx    = grid spacing in x
% dz    = grid spacing in z
% phi   = nx by nz matrix containing values of phi
% phi_x = nx by nz matrix containing gradient of phi in x-direction at each cell center
% phi_z = nx by nz matrix containing gradient of phi in z-direction at each cell center

phi_x = zeros(nx,nz);
phi_z = zeros(nx,nz);

for i = 1:nx
    for j = 1:nz
        
        % x-gradient:
        if i==1
            phip = phi(i+1,j);
            phim = phi(nx,j);        %based on periodic condition with FVM
        elseif i==nx
            phip = phi(1,j);
            phim = phi(i-1,j);
        else
            phip = phi(i+1,j);
            phim = phi(i-1,j);
        end
        phi_x(i,j) = (phip-phim)/dx/2.0;
        
        
        % z-gradient:
        if j==1
            phip = phi(i,j+1);
            phim = phi(i,nz);
        elseif j==nz
            phip = phi(i,1);
            phim = phi(i,j-1);
        else
            phip = phi(i,j+1);
            phim = phi(i,j-1);
        end
        phi_z(i,j) = (phip-phim)/dz/2.0;
        
    end
end


end

