clear, clc

%Variables%
Theta_l = ; % Motor angle left
Theta_r = ; % Motor angle right
u_ml = ;    % Motor voltage left
u_mr = ;    % Motor voltage right
i_ml = ;    % Motor current left
i_mr = ;    % Motor current right
r_m = ;     % Motor radius
L_m = ;     % Motor induction
R_m = ;     % Motor resistor
J_m = ;     % Motor inertia
K_m = ;     % Motor constant
m_tool = ;  % Tool head mass
u_n = ;     % Nozzle voltage
i_n = ;     % Nozzle current
C_n = ;     % Nozzle capacitor
C_tn = ;    % Nozzle thermal capacity
R_tn = ;    % Nozzle thermal resistance
R_n = ;     % Nozzle resistor
T_n = ;     % Nozzle temperature
T_e = ;     % Environmental temperature
MU_r = ;    % Bearing friction rolling
MU_lx = ;   % Bearing friction linear x-direction
MU_ly = ;   % Bearing friction linear y-direction
%Theta_l = ;     % Motor angle left
%Theta_r = ;     % Motor angle right
%u_ml = ;        % Motor voltage left
%u_mr = ;        % Motor voltage right
%i_ml = ;        % Motor current left
%i_mr = ;        % Motor current right
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
%T_n = ;         % Nozzle temperature
%T_e = ;         % Environmental temperature
mu_r = 0.3;        % Bearing friction rolling [N*s/m]
mu_lx = 0.5;       % Bearing friction linear x-direction [N*s/m]
mu_ly = 0.5;       % Bearing friction linear y-direction [N*s/m]
%T_b = ;         % Print bed temperature
dx = 0.3;         % Print area dimension x [m]
dy = 0.3;         % Print area dimension y [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_eff = 1/(J_m + 0.5*m_tool*r_m^2);   % effective inertia

% x-dynamics
A_x = [0,               1,                                  0;
       0,  -J_eff * r_m^2 * (4*mu_r + mu_lx),     (J_eff*r_m*K_m)/2;
       0,       -(2*K_m)/(r_m*L_m),                    -R_m/L_m];
B_x = [0; 0; 1/L_m];
C_x = [1, 0, 0];
D_x = 0;

sys_x = ss(A_x, B_x, C_x, D_x);

% y-dynamics
A_y = [0,               1,                                  0;
       0,  -J_eff * r_m^2 * (2*mu_r + mu_ly),     (J_eff*r_m*K_m)/2;
       0,       -(2*K_m)/(r_m*L_m),                    -R_m/L_m];
B_y = [0; 0; 1/L_m];
C_y = [1, 0, 0];
D_y = 0;

sys_y = ss(A_y, B_y, C_y, D_y);

% Inputs
t    = 0:0.001:5;                % Simulation time [s]
f    = 0.125;                     % Frequency [Hz] (1 circle per second)
w    = 2 * pi * f;            % Angular frequency [rad/s]
U_0  = 1.5;                    % Input voltage amplitude [V]

% Decoupled inputs
u_x = U_0 * cos(2* w * t);
%u_x = 0*t;
u_y = U_0 * sin(2* w * t);
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

% Calculate actual physical motor voltages
u_ml = (u_x + u_y) / 2;
u_mr = (u_x - u_y) / 2;

% Run 
x_pos = lsim(sys_x, u_x, t);
y_pos = lsim(sys_y, u_y, t);

% Initialize combined figure
figure('Name', 'CoreXY Toolhead Motion', 'Color', 'w');

% Position vs. Time
subplot(1, 2, 1);
plot(t, x_pos, 'b-', 'LineWidth', 1.5, 'DisplayName', 'x(t)');
hold on;
plot(t, y_pos, 'r--', 'LineWidth', 1.5, 'DisplayName', 'y(t)');
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
xlabel('X Position [m]');
ylabel('Y Position [m]');
title('Physical Toolpath Trajectory');