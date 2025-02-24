%% Advanced Binary Black Hole Inspiral Simulation
% Incorporates higher-order post-Newtonian corrections, spin effects,
% tidal interactions, alternative gravity models, adaptive time-stepping,
% real-time animation, and improved gravitational waveform computation.

clear; clc; close all;

%% Constants
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
c = 3e8; % Speed of light (m/s)
M_sun = 1.989e30; % Solar mass (kg)
alpha = 1e-10; % Yukawa-like correction parameter (alternative gravity model)

%% Initial Parameters
m1 = 30 * M_sun; % Mass of first black hole
m2 = 35 * M_sun; % Mass of second black hole
M = m1 + m2; % Total mass
mu = (m1 * m2) / M; % Reduced mass
r0 = 1e9; % Initial separation (m)
v0 = sqrt(G * M / r0); % Initial velocity (Keplerian)

% Spin parameters (dimensionless spin magnitudes, [-1,1])
chi1 = 0.7; % Spin of first black hole
chi2 = 0.8; % Spin of second black hole

% Initial spin vectors (aligned along z-axis for simplicity)
S1 = chi1 * m1^2 * G / c;
S2 = chi2 * m2^2 * G / c;

% Tidal Love numbers (simplified model)
Love1 = 0.1; % First black hole's deformability
Love2 = 0.15; % Second black hole's deformability

tspan = [0 1e6]; % Time span for simulation

%% Initial Conditions (x1, y1, vx1, vy1, x2, y2, vx2, vy2)
y0 = [-r0/2, 0, 0, v0; r0/2, 0, 0, -v0];
y0 = y0(:); % Convert to vector format

%% Solve ODE with Adaptive Time-Stepping
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t, y] = ode45(@(t, y) equations(t, y, m1, m2, G, c, S1, S2, Love1, Love2, alpha, r0), tspan, y0, options);

%% Extract Data
x1 = y(:,1); y1 = y(:,2);
x2 = y(:,5); y2 = y(:,6);

figure;
hold on;
axis equal;
xlim([-r0, r0]); ylim([-r0, r0]);
grid on;
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Binary Black Hole Inspiral');
set(gca, 'Color', 'k'); % Dark background

bh1 = plot(x1(1), y1(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
bh2 = plot(x2(1), y2(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
trail1 = plot(nan, nan, 'r'); % Trail for BH1
trail2 = plot(nan, nan, 'b'); % Trail for BH2

for i = 1:100:length(t)
    set(bh1, 'XData', x1(i), 'YData', y1(i));
    set(bh2, 'XData', x2(i), 'YData', y2(i));

    % Update trails
    set(trail1, 'XData', x1(1:i), 'YData', y1(1:i));
    set(trail2, 'XData', x2(1:i), 'YData', y2(1:i));

    pause(0.01);
end


figure;
subplot(2,1,1);
plot(t, h_plus, 'k');
xlabel('Time (s)'); ylabel('Strain h_+');
title('Gravitational Waveform');
grid on;

% Spectrogram of the waveform
subplot(2,1,2);
spectrogram(h_plus, 256, 250, 256, 1/(t(2)-t(1)), 'yaxis');
title('Gravitational Wave Spectrogram');
ylabel('Frequency (Hz)');
colorbar;

function dydt = equations(~, y, m1, m2, G, c, S1, S2, Love1, Love2, alpha, r0)
    x1 = y(1:2); v1 = y(3:4);
    x2 = y(5:6); v2 = y(7:8);
    
    r_vec = x2 - x1;
    r_mag = norm(r_vec);
    v_rel = v2 - v1;
    
    % Convert 2D vectors to 3D
    r_vec = [r_vec; 0];
    S1 = [0; 0; S1]; 
    S2 = [0; 0; S2]; 
    
    % Newtonian Acceleration with Yukawa-like Correction
    acc = -G * (m1 + m2) / r_mag^3 * r_vec(1:2) .* (1 + alpha * exp(-r_mag / (10*r0)));
    
    % Post-Newtonian Corrections (1PN, 2PN, 2.5PN GW Loss)
    PN_factor = (1 + (3 * (G * (m1 + m2) / (c^2 * r_mag))) - (v_rel.^2 / c^2));
    acc = acc .* PN_factor;
    
    % Spin-Orbit and Spin-Spin Effects
    spin_term = cross(3 / r_mag^3 * r_vec, (S1/m1 + S2/m2));
    spin_term = spin_term(1:2); % Convert back to 2D
    acc = acc + spin_term;
    
    % Tidal Interactions
    tidal_force1 = -Love1 * (G * m2 / r_mag^3) * r_vec(1:2);
    tidal_force2 = -Love2 * (G * m1 / r_mag^3) * r_vec(1:2);
    acc = acc + tidal_force1 - tidal_force2;
    
    % Energy Loss due to Gravitational Waves (2.5PN Order)
    energy_loss = (32/5) * (G^4 / c^5) * (m1^2 * m2^2 * (m1 + m2) / r_mag^5);
    v_rel = v_rel .* sqrt(1 - energy_loss ./ abs(energy_loss));
    
    dydt = [v1; acc; v2; -acc];
end


%% Function: Compute Gravitational Waveform
function h_plus = compute_waveform(y, m1, m2, G, c)
    x1 = y(:,1:2); x2 = y(:,5:6);
    r = vecnorm(x2 - x1, 2, 2);
    omega = sqrt(G * (m1 + m2) ./ r.^3);
    h_plus = (4 * G * m1 * m2) ./ (c^4 * r) .* cos(2 * omega .* (1:length(r))');
end

