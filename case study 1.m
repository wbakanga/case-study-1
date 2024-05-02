% Parameters
R_inner = 0.1;      % Inner radius of the pipe (m)
R_outer = 0.15;     % Outer radius of the pipe (m)
L = 5;              % Length of the pipe (m) - increased length
Nr = 50;            % Number of radial grid points
Ntheta = 50;        % Number of angular grid points
dr = (R_outer - R_inner) / (Nr - 1); % Radial grid spacing
dtheta = 2 * pi / (Ntheta - 1);      % Angular grid spacing
alpha = 0.01;       % Thermal diffusivity (m^2/s)
T_inner = 100;      % Inner surface temperature (°C) - increased temperature
T_outer = 20;       % Outer surface temperature (°C)
dt = 0.001;         % Time step (s)
t_final = 1;        % Final simulation time (s)

% Initial conditions
T0 = T_outer * ones(Nr, Ntheta); % Initial temperature

% Main loop
for t = 0:dt:t_final
    T_new = T0;
    for i = 2:Nr-1
        for j = 2:Ntheta-1
            % Finite difference scheme for heat conduction
            T_new(i, j) = T0(i, j) + alpha * dt * ((T0(i+1, j) - 2*T0(i, j) + T0(i-1, j)) / dr^2 + ...
                (1/i) * (T0(i+1, j) - T0(i-1, j)) / (2 * dr) + (1/(i^2)) * (T0(i+1, j) - 2*T0(i, j) + T0(i-1, j)) / dr^2 + ...
                (T0(i, j+1) - 2*T0(i, j) + T0(i, j-1)) / (i * dtheta * dr));
        end
    end
    
    % Boundary conditions (Dirichlet)
    T_new(1, :) = T_inner;   % Inner surface temperature
    T_new(Nr, :) = T_outer;  % Outer surface temperature
    
    T0 = T_new; % Update temperature matrix
end

% Plotting 2D contour
theta = linspace(0, 2*pi, Ntheta);
r = linspace(R_inner, R_outer, Nr);
[R, Theta] = meshgrid(r, theta);
X = R .* cos(Theta);
Y = R .* sin(Theta);
contourf(X, Y, T0', 'LineStyle', 'none');
colormap(jet);
colorbar;
xlabel('x');
ylabel('y');
title('Heat Conduction through Pipe from Internal to External Surface (Contour Plot)');