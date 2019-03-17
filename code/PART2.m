% ELEC 4700 Assignment 3 Part 2
% Andrew Branicki 100973961
% March 17 2019

global C
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons
    
% Declare boundaries of the canvas
% For plots, use a ratio of 3/2 for L/W
% Similar to EIPA from February 6, 2019
L = 150; W = 100; dx = 1; dy = 1;
G = sparse(L*W, W*L);
V = zeros(L*W,1);
cond = zeros(L,W);
map = @(i,j) j + (i - 1)*W;

% Conductivity values
sigma_o = 1;        % outside boxes 
sigma_in = 1e-2;    % inside boxes

% set up the regions for the boxes
region_x = L/3;
region_y = W/3;
middle_x = L/2;
middle_y = W/2;

% Set up the conductivity map
for i = 1:L
    for j = 1:W
        % check if we are within the boxes
        if ((i > (middle_x - region_x/2) && i < (middle_x + region_x/2)) && (j > (middle_y + region_y/2) || j < (middle_y - region_y/2)))  
            cond(i,j) = sigma_in;
        % else we are not within the boxes
        else
            cond(i,j) = sigma_o;
        end
    end
end

% X direction
for i=1:L
    % Y direction
    for j=1:W
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        % nx = 0 (in this case, V = Vo)
        % we also don't have boxes on the left and right edges so don't
        % need to check anything :)
        if i == 1
            G(n,:) = 0;
            G(n,n) = cond(i,j);
            V(n) = 1;
        % nx = L (in this case, V = 0)
        elseif i == L
            G(n,:) = 0;
            G(n,n) = cond(i,j);
            V(n) = 0;
        % Y edges of canvas
        elseif j == 1
            G(n,:) = 0;
            G(n,nxm) = (cond(i-1,j) + cond(i,j))/2;
            G(n,nxp) = (cond(i+1,j) + cond(i,j))/2;
            G(n,nyp) = (cond(i,j+1) + cond(i,j))/2;
            G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nyp));
        elseif j == W
            G(n,:) = 0;
            G(n,nxm) = (cond(i-1,j) + cond(i,j))/2;
            G(n,nxp) = (cond(i+1,j) + cond(i,j))/2;
            G(n,nym) = (cond(i,j-1) + cond(i,j))/2;
            G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nym));
        
        else          
            G(n,:) = 0;
            G(n,nxm) = (cond(i-1,j) + cond(i,j))/2;
            G(n,nxp) = (cond(i+1,j) + cond(i,j))/2;
            G(n,nyp) = (cond(i,j+1) + cond(i,j))/2;
            G(n,nym) = (cond(i,j-1) + cond(i,j))/2;
            G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nyp) + G(n,nym));  
        end

    end
end

% Solve for F using G\V = F
F = G\V;

% Set up a surf plot
surfs_up = zeros(L,W);
for i = 1:L
    for j = 1:W
        n = map(i,j);
        surfs_up(i,j) = F(n);
    end
end

% TIME TO PLOT COOL THINGS!
[dimensions_X, dimensions_Y] = meshgrid(1:W,1:L);
[EX, EY] = gradient(surfs_up);

figure(5)
surf(surfs_up)
colorbar
title('Voltage Map')
xlabel('Height (nm)')
ylabel('Width(nm)')
zlabel('Voltage (V)')
view(135,45)

figure(6)
quiver(dimensions_X, dimensions_Y, EX, EY)
title('Electric Field Quiver Plot')
xlabel('Height (nm)')
ylabel('Width(nm)')
