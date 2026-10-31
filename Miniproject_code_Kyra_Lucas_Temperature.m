clear, clc

%Variables%
r_m = 0.018;         % Motor radius [m]
L_m = 1e-4;         % Motor induction [H]
R_m = 2.0;         % Motor resistor [Ohm]
J_m = 1e-5;         % Motor inertia [kg*m^2]
K_m = 0.05;         % Motor constant [N*m/A]
m_tool = 0.5;      % Tool head mass [kg]
%u_n = ;         % Nozzle voltage
%i_n = ;         % Nozzle current
%C_n = ;         % Nozzle capacitor
%C_tn = ;        % Nozzle thermal capacity
%R_tn = ;        % Nozzle thermal resistance
%R_n = ;         % Nozzle resistor
mu_r = 0.3;        % Bearing friction rolling [N*s/m]
mu_lx = 0.5;       % Bearing friction linear x-direction [N*s/m]
mu_ly = 0.5;       % Bearing friction linear y-direction [N*s/m]
dx = 0.3;         % Print area dimension x [m]
dy = 0.3;         % Print area dimension y [m]

% Heat Transfer & Thermal Variables %
d = 0.0004;        % Filament cross-section dimension [m]
rho = 1240;        % Filament density [kg/m^3]
c_p = 1800;        % Specific heat capacity [J/(kg*K)]
k_cond = 0.13;     % Thermal conductivity [W/(m*K)]
alpha = k_cond/(rho*c_p); % Thermal diffusivity [m^2/s]
h_air = 20;        % Air convection coefficient [W/(m^2*K)]
h_bed = 100;       % Bed contact conduction coefficient [W/(m^2*K)]
T_e = 25;          % Environmental temperature [°C]
T_b = 60;          % Print bed temperature [°C]
T_n = 210;         % Nozzle temperature [°C]
N = 100;            % Number of spatial grid cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_eff = 1/(J_m + 0.5*m_tool*r_m^2);   % effective inertia

% x-dynamics
A_x = [0,               1,                                  0;
       0,  -J_eff * r_m^2 * (4*mu_r + mu_lx),     (J_eff*r_m*K_m)/2;
       0,       -(2*K_m)/(r_m*L_m),                    -R_m/L_m];
B_x = [0; 0; 1/L_m];
C_x = [1, 0, 0];
D_x = 0;

eig(A_x)
sys_x = ss(A_x, B_x, C_x, D_x);


% Inputs
t    = 0:0.001:15;                % Simulation time [s]
f    = 0.025;                     % Frequency [Hz] (1 circle per second)
w    = 2 * pi * f;            % Angular frequency [rad/s]
U_0  = 0.05;                    % Input voltage amplitude [V]

% Decoupled inputs
u_x = U_0 * cos(5* w * t);
%u_x = 0*t;
u_y = U_0 * sin(5* w * t);
%u_y = 1*t;

if u_x == 0
    A_y = [0,               1,                                  0;
           0,  -J_eff * r_m^2 * (2*mu_r + mu_ly),     (J_eff*r_m*K_m)/2;
           0,       -(2*K_m)/(r_m*L_m),                    -R_m/L_m];
    B_y = [0; 0; 1/L_m];
    C_y = [1, 0, 0];
    D_y = 0;

    sys_y = ss(A_y, B_y, C_y, D_y);

else
    A_y = [0,               1,                                  0;
       0,  -J_eff * r_m^2 * (4*mu_r + mu_ly),     (J_eff*r_m*K_m)/2;
       0,       -(2*K_m)/(r_m*L_m),                    -R_m/L_m];
    B_y = [0; 0; 1/L_m];
    C_y = [1, 0, 0];
    D_y = 0;

    sys_y = ss(A_y, B_y, C_y, D_y);
end

eig(A_y)
% Calculate actual physical motor voltages
u_ml = (u_x + u_y) / 2;
u_mr = (u_x - u_y) / 2;

% Run 
x_pos = lsim(sys_x, u_x, t);
y_pos = lsim(sys_y, u_y, t);

% Initialize combined figure
figure('Name', 'CoreXY Toolhead Motion');

% Position vs. Time
subplot(1, 2, 1);
plot(t, x_pos, 'b-', 'DisplayName', 'x(t)');
hold on;
plot(t, y_pos, 'r--', 'DisplayName', 'y(t)');
hold off;
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
title('Positions vs. Time');
legend('Location', 'best');

% 2D Print Bed Trajectory (XY Plane)
subplot(1, 2, 2);
plot(x_pos, y_pos, 'm-', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('X Position [m]');
ylabel('Y Position [m]');
title('Tool head position');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$ Q1.6

%% Calculate Path Length and Print Speed
T_loop = 1 / (5 * f);
N_loop = round(T_loop / 0.001);

x1 = x_pos(1:N_loop);
y1 = y_pos(1:N_loop);

% Distance between consecutive path points
dist = sqrt(diff(x1).^2 + diff(y1).^2);
L = sum(dist);              % Total loop length [m]
v = mean(dist / 0.001);     % Average print speed [m/s]
x_grid = linspace(0, L*1000, N);

% Thermal Model Setup
dl = L / N;                 % Distance between grid nodes
t_print = L / v;            % Print duration [s]
t_final = t_print + 25;     % Total simulation time (print + cooling) [s]

dt = (dl / v) / 10;         % Stable time step [s]
k_f = floor(t_final / dt) + 1;

M = NaN(N, k_f);          % Initialize history matrix with NaN
T = NaN(N, 1);            % Unprinted filament is NaN (does not exist yet)
T(1) = T_n;               % First extruded point
M(:, 1) = T;

n_i = 1;                    % Currently extruded node count

% Thermal Simulation Loop
for k = 1:(k_f - 1)
    t_sim = (k - 1) * dt;
    
    if t_sim <= t_print
        % PHASE 1: Extruding filament (growing domain)
        n_i_new = floor(v * t_sim / dl) + 1;
        if n_i_new > n_i && n_i_new <= N
            T(n_i:n_i_new) = T_n;
            n_i = n_i_new;
        end
        
        T_next = T;
        for i = 2:(n_i - 1)
            dT = alpha/(dl^2) * (T(i-1) + T(i+1) - 2*T(i)) + 3*h_air/(rho*c_p*d)*(T_e - T(i)) + h_bed/(rho*c_p*d)*(T_b - T(i));
            T_next(i) = T(i) + dT * dt;
        end
        T = T_next;

        if n_i >= 2
            T(1) = T(2);
        end
        
    else
        % PHASE 2: Closed loop cooling (hot end connects to cold start)
        for i = 1:N
            if i == 1
                T_left = T(N);   T_right = T(2);
            elseif i == N
                T_left = T(N-1); T_right = T(1);
            else
                T_left = T(i-1); T_right = T(i+1);
            end
            
            dT = alpha/(dl^2) * (T_left + T_right - 2*T(i)) + 3*h_air/(rho*c_p*d)*(T_e - T(i)) + h_bed/(rho*c_p*d)*(T_b - T(i));
            T_next(i) = T(i) + dT * dt;
        end
        T = T_next;
    end
    
    M(:, k+1) = T; % Store current temperature column
end

% Plots
t_vec = (0:(k_f-1)) * dt;       % time array for heatmap [s]

figure('Name', 'Filament Temperature Analysis', 'Position', [100, 100, 1000, 450]);

% Subplot 1: Temperature Profile Animation along filament
subplot(1, 2, 1);
for k = 1:20:k_f
    plot(x_grid, M(:, k), 'b-', 'LineWidth', 1.8);
    title(['Time ' num2str((k-1)*dt, '%.2f') ' [s]']);
    xlabel('Position along filament [mm]');
    ylabel('Temperature [°C]');
    ylim([20 220]);
    grid on;
    pause(1/120);
end

% Subplot 2: Spatiotemporal Heat Map (Time vs. Position)
subplot(1, 2, 2);
imagesc(t_vec, x_grid, M);
set(gca, 'YDir', 'normal');          % Keep position 0 at bottom
colormap(jet);
c = colorbar;
c.Label.String = 'Temperature [°C]';
xlabel('Time [s]');
ylabel('Position along filament [mm]');
title('Filament Heat Map');
xline(t_print, 'w--', 'Print Finished', 'LineWidth', 1.5, 'LabelOrientation', 'aligned');

% Find the time column index at the moment printing finishes
k_print = floor(t_print / dt) + 20;

% Plot temperature along the filament at t = t_print
figure('Name', 'Temperature at Print Completion');
plot(x_grid, M(:, k_print), 'b-', 'LineWidth', 1.8);
grid on;
xlabel('Position along filament [mm]');
ylabel('Temperature [°C]');
title(['Temperature Profile at Print Completion (t = ' num2str(t_print, '%.2f') ' s)']);
ylim([20 220]);
