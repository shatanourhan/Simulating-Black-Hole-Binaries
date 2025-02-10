clc; clear; close all;

% Constants
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
c = 3.0e8;       % Speed of light (m/s)
m1 = 10e30;      % Mass of black hole 1 (kg)
m2 = 10e30;      % Mass of black hole 2 (kg)
M = m1 + m2;     % Total mass

% Initial conditions
r0 = [1e9; 0];   % Initial position (m)
v0 = [0; 1e4];   % Initial velocity (m/s)
y0 = [r0; v0];

tspan = [0 1e6]; % Time span (s)
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Solving the system with post-Newtonian corrections
theta = linspace(0, 2*pi, 500);
[t, y] = ode45(@(t, y) equations(t, y, G, m1, m2, c), tspan, y0, options);

% Plot results
figure;
plot(y(:,1), y(:,2), 'b'); hold on;
plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Center of mass
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Binary Black Hole Inspiral with Post-Newtonian Effects');
grid on;
legend('Orbit', 'Center of Mass');

% Computing and plot gravitational waveform
h_plus = compute_waveform(y, G, m1, m2, c);
figure;
plot(t, h_plus);
xlabel('Time (s)'); ylabel('Strain h_+');
title('Gravitational Waveform');
grid on;

function dydt = equations(~, y, G, m1, m2, c)
    r = y(1:2);
    v = y(3:4);
    r_mag = norm(r);
    
    % Newtonian acceleration
    acc = -G * (m1 + m2) * r / r_mag^3;
    
    % Post-Newtonian corrections (1PN, 2PN, 2.5PN for GW loss)
    v_mag = norm(v);
    PN_factor = 1 + (3 * G * (m1 + m2) / (c^2 * r_mag)) - (v_mag^2 / c^2);
    acc = acc * PN_factor;
    
    % Gravitational wave energy loss (2.5PN order)
    power = (32/5) * (G^4 / c^5) * (m1^2 * m2^2 * (m1 + m2) / r_mag^5);
    energy = 0.5 * (m1 + m2) * norm(v)^2 - G * m1 * m2 / r_mag;
    
    if energy < 0
        v = v * sqrt(1 - power / abs(energy));
    end
    
    dydt = [v; acc];
end

function h_plus = compute_waveform(y, G, m1, m2, c)
    r = vecnorm(y(:,1:2), 2, 2);
    omega = sqrt(G * (m1 + m2) ./ r.^3);
    
    % Quadrupole formula for GW strain
    h_plus = (4 * G * m1 * m2) ./ (c^4 * r) .* cos(2 * omega .* (1:length(r))');
end
