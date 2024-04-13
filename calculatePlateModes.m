function [fLong, fTrans,fBend,fEig] = calculatePlateModes(length, width, thickness, n_modes, density, E, P)
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
    %   fLong     - Longitudinal modes of vibration for the plate as a 2 by
    %               n_modes vector representing the longitudinal modes for
    %               length and width
    %   fTrans    - Transverse shear modes of vibration for the plate as a 2 
    %               by n_modes vector representing the transverse modes for
    %               length and width
    %   fBend     - Bending modes of vibration for the plate as a 2 by n_modes
    %               vector representing the flexural modes for length and width
    %   fEig      - The modes found from the eigenvalues generated from the 
    %               finite difference method
    
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
    % [longitudinal_modes,transversal_modes,beinding_modes,eigenvalue_modes] = calculatePlateModes(length ...
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
    % Define the plate modes by replicating bending modes
    fPlate = zeros(2, n_modes);
    for i = 1:n_modes
        fPlate(1, i) = fBend(1, i); % Along length
        fPlate(2, i) = fBend(2, i); % Along width
    end
   
    
    %% Modal Analysis (Extrapolation to Plate Modes)
    % Define the plate modes by replicating bending modes
    fPlate = zeros(2, n_modes);
    for i = 1:n_modes
        fPlate(1, i) = fBend(1, i); % Along length
        fPlate(2, i) = fBend(2, i); % Along width
    end
    
    %% Finite Element Analysis
    
    % Define number of grid points and grid spacing
    Nx = 30; 
    Ny = 30; 
    h = length / (Nx - 1); 
    k = width / (Ny - 1); 
    
    % Create Grid
    x = linspace(0, length, Nx);
    y = linspace(0, width, Ny);
    [X, Y] = meshgrid(x, y);
    
    % Formulate the stiffness matrix (using finite differences approximation)
    K = zeros(Nx*Ny, Nx*Ny);
    
    for i = 1:Nx
        for j = 1:Ny
            % Node number
            n = (j - 1) * Nx + i;
            
            % Calculate plate equation terms
            Dxx = E * thickness^3 / (12 * (1 - P^2));
            Dyy = E * thickness^3 / (12 * (1 - P^2)); 
            
            % Define terms for stiffness matrix
            if i > 1
                % Left neighbor
                K(n, n-1) = Dxx / (h^4);
            end
            if i < Nx
                % Right neighbor
                K(n, n+1) = Dxx / (h^4);
            end
            if j > 1
                % Bottom neighbor
                K(n, n-Nx) = Dyy / (k^4);
            end
            if j < Ny
                % Top neighbor
                K(n, n+Nx) = Dyy / (k^4);
            end
            % Self term
            K(n, n) = -2 * (Dxx / (h^4) + Dyy / (k^4));
        end
    end

    
    % Solve the eigenvalue problem using sparse solver
    num_modes_to_compute = n_modes;
    % Suppress output
    options = struct('disp',0); 
    
    % Initialize variables for convergence check
    
        % Solve eigenvalue problem
    [V, D] = eigs(K, num_modes_to_compute, 'SM', options);
    %[V, D] = eig(K);
        
    % Extract eigenvalues (frequencies) and eigenvectors (mode shapes)
    % Convert angular frequencies to Hz
    fEig = diag(sqrt(-D./(2*pi))); 
    
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
        title(['Mode ', num2str(i)]);
        colorbar;
    end
    end
    