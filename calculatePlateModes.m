function [fLong, fTrans,fBend] = calculatePlateModes(length, width, thickness, n_modes, density, E, P)
%   This function is meant to give different approaches to the calculation
%   of modes for a plate.there is a computational method, which is
%   the returned frequencies, and there is the eigenvalue method which is 
%   generated and visualized according to classical plate theory.
%
% Inputs:
%   length    - Length of the steel plate in m
%   width     - Width of the steel plate in m
%   thickness - Thickness of the steel plate in m
%   n_modes   - Number of modes you want to calculate
%   density   - Density of the material in kg m^-3
%   E         - Youngs modulus of the material in GPa
%   P         - Poisson ratio of the material
%
% Outputs:
%   fLong     - Longitudinal modes of vibration for the plate
%   fTrans    - Transverse shear modes of vibration for the plate
%   fBend     - Bending modes of vibration for the plate

%% Example for plate of 4x8 feet made from 26 gauge galvanized steel

%length = 1.2192;
% width = 2.4384;
% thickness = 0.00045;
% n_modes = 100;
% density = 7800;
% youngs_mod = 210;
% poisson = 0.29;
% 
% 
% [longitudinal_modes,transversal_modes,beinding_modes] = calculatePlateModes(length ...
%     ,width,thickness,n_modes,density,youngs_mod,poisson);


%% Start of code
% unit conversion
E = E * 1e9;

%% Longitudinal
% defining speed of sound of longitudinal waves in plate
cLong = sqrt(E / density);
% initializing frequency vector
fLong = zeros(2, n_modes);
for n = 1:n_modes
    fLong(1, n) = n * cLong / length;
    fLong(1, n) = round(fLong(1, n));
    fLong(2, n) = n * cLong / width;
    fLong(2, n) = round(fLong(2, n));
end

%% Transverse shear
% defining speed of sound of transverse shear waves in plate
G = E / (2 * (1 + P));
cTrans = sqrt(G / density);
fTrans = zeros(2, n_modes);
for n = 1:n_modes
    fTrans(1, n) = n * cTrans / (length ^ 2);
    fTrans(1, n) = round(fTrans(1, n));
    fTrans(2, n) = n * cTrans / (width ^ 2);
    fTrans(2, n) = round(fTrans(2, n));
end


%% Bending modes
fBend = zeros(2, n_modes);
for n = 1:n_modes
    fBend(1, n) = (n^2 * pi^2 * E * thickness / (12 * (1 - P^2) * density))^0.5 / length^2;
    fBend(1, n) = round(fBend(1, n));
    fBend(2, n) = (n^2 * pi^2 * E * thickness / (12 * (1 - P^2) * density))^0.5 / width^2;
    fBend(2, n) = round(fBend(2, n));
end


%% Classical plate theory

% Define number of grid points and grid spacing
Nx = 30; 
Ny = 30; 
h = length / (Nx - 1); 
k = width / (Ny - 1); 

% Initialize deflection matrix
w = zeros(Ny, Nx);

% Create Grid
x = linspace(0, length, Nx);
y = linspace(0, width, Ny);
[X, Y] = meshgrid(x, y);

% Apply Boundary Conditions (Simply Supported)
w(:, 1) = 0;
w(:, end) = 0;
w(1, :) = 0;
w(end, :) = 0;

% Solve Plate Equation Using Finite Differences (Without Load)
for iter = 1:10000 %(adjust or use convergence criteria)
    w_old = w; % Store old deflection for convergence check
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Plate equation in 2D
            wxx = (w(j, i+1) - 2*w(j, i) + w(j, i-1)) / h^2;
            wyy = (w(j+1, i) - 2*w(j, i) + w(j-1, i)) / k^2;
            w(j, i) = (E*thickness^3 / (12*(1 - thickness^2))) * (wxx + 2*(1 + thickness)*wyy);
        end
    end
    
    % Check for convergence
    if max(abs(w_old(:) - w(:))) < 1e-6
        break;
    end
end

%% Modal Analysis
% Formulate the stiffness matrix (using finite differences approximation)
K = zeros(Nx*Ny);
for i = 1:Nx
    for j = 1:Ny
        % Node number
        n = (j - 1) * Nx + i; 
        if i > 1
            % Left neighbor
            K(n, n-1) = 1 / h^2; 
        end
        if i < Nx
            % Right neighbor
            K(n, n+1) = 1 / h^2; 
        end
        if j > 1
            %Bottom neighbor
            K(n, n-Nx) = 1 / k^2;
        end
        if j < Ny
            % Top neighbor
            K(n, n+Nx) = 1 / k^2; 
        end
        % Self term
        K(n, n) = -2 * (1/h^2 + 1/k^2); 
    end
end

% Solve the eigenvalue problem using sparse solver
num_modes_to_compute = 9;
% Suppress output
options = struct('disp',0); 
% 'SM' for smallest magnitude eigenvalues
[V, D] = eigs(K, num_modes_to_compute, 'SM', options); 

% Extract eigenvalues (frequencies) and eigenvectors (mode shapes)
% Convert angular frequencies to Hz
frequencies = sqrt(-D)./2*pi; 

% Plot Mode Shapes
num_modes_to_plot = 9;
figure;
for i = 1:num_modes_to_plot
    mode_shape = reshape(V(:, i), Ny, Nx);
    subplot(3, 3, i);
    surf(X, Y, mode_shape);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('Deflection (m)');
    title(['Mode at ', num2str(frequencies(i,i)), ' Hz']);
    colorbar;
end

end
