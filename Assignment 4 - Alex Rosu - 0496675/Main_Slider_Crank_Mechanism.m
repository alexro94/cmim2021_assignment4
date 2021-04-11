%% Assignment 4 - Main Program
% Alex Rosu
% 11.04.2021

clear; clc;

% u = [theta; d];
a = 0.1;                    % The length of the crank [m]
b = 0.2;                    % The length of the shaft [m]
phi = pi/6;                 % Initial angular position of the crank [rad]
omega = 1;                  % Angular velocity of the crank [rad/s]
phi_dot = -1;               % Angular velocity opposite to omega [rad/s]

t = 0:.1:10;                % Simulation time
pos = zeros(length(t), 2);  % Initializing vector for collecting positions
vel = zeros(length(t), 2);  % Initializing vector for collecting velocities

eps = 1e-9;                 % allowed difference between exact value and approx.

% set a reasonable starting point for positions and velocities
u0 = [0.25; b + a];
u0_der = [-0.4; 0.07];

% Determining positions and velocities in a for-loop
for i = 1:length(t)
    
    phi = pi/6 + omega*t(i);   % Updating angle phi
    
    % Create function handles to determine theta and d
    F = @(u) constraint(u, a, b, phi);
    J = @(u) jacobian(u, b);
    
    % Calling Newton-Raphson function to determine valid position
    [u, iteration_counter] = NR_method(F, J, u0, eps);
    pos(i,:) = u; % Updating position vector
    
    % Create function handles to determine theta_dot and d_dot
    F_der = @(u_der) constraint_der(u_der, a, b, phi, phi_dot, u(1));
    J_der = @(u_der) jacobian_der(b, u(1));

    % Calling Newton-Raphson function to determine valid velocities
    [u_der, iteration_counter_der] = NR_method(F_der, J_der, u0_der, eps);
    vel(i,:) = u_der; % Updating velocity vector
     
end

% Plotting the results
plot(t, pos(:,1), 'r');
grid on;
xlabel('Time [s]');
ylabel('Angular position of the shaft [rad]');
axis tight;
set(gca,'FontSize',14)

figure(2);
plot(t, pos(:,2), 'r');
grid on;
xlabel('Time [s]');
ylabel('Position of the block [m]');
axis tight;
set(gca,'FontSize',14)

figure(3);
plot(t, vel(:,1), 'r');
grid on;
xlabel('Time [s]');
ylabel('Angular velocity of the shaft [rad/s]');
axis tight;
set(gca,'FontSize',14)

figure(4);
plot(t, vel(:,2), 'r');
grid on;
xlabel('Time [s]');
ylabel('Velocity of the block [m/s]');
axis tight;
set(gca,'FontSize',14)


% Equations for position constraints
function P = constraint(u, a, b, phi)
theta = u(1);
d = u(2);
P = [a * cos(phi) + b * cos(theta) - d
    a * sin(phi) - b * sin(theta)];
end

% Jacobian matrix for position constraints
function P = jacobian(u, b)
theta = u(1);
P = [-b * sin(theta), -1
    -b * cos(theta), 0];
end

% Equations for velocity constraints
function P_der = constraint_der(u_der, a, b, phi, phi_dot, theta)
theta_dot = u_der(1);
d_dot = u_der(2);
P_der = [-a * phi_dot * sin(phi) - b * theta_dot * sin(theta) - d_dot
        a * phi_dot * cos(phi) - b * theta_dot * cos(theta)];
end

% Jacobian matrix for velocity constraints
function P_der = jacobian_der(b, theta)
P_der = [-b * sin(theta), -1
        -b * cos(theta), 0];
end