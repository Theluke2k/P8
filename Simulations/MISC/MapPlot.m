%% Initial Conditions Map Plot
% This script visualizes the initial map: process plume, robots, and charging station

% Grid for map
xmin = 0; xmax = 40;
ymin = 0; ymax = 40;
zmin = 0; zmax = 5;
x_axis = linspace(xmin, xmax, 100);
y_axis = linspace(ymin, ymax, 100)';

% Initial process state (scaled back to real units)
sc = [200 1 0.05 0.0001 20 0.1 20 0.1]';
T = diag(sc);
z0 = [300; 0; 0.05; 0; 20; 0; 20; 0];       % true initial vector
z0_scaled = inv(T)*z0;

% Compute process height at grid
[Xg, Yg] = meshgrid(x_axis, y_axis);
Z_true = get_h(z0_scaled, x_axis, y_axis, sc);

% Charging station parameters
t1 = 0.5; % sharpness
h = 5;    % half-size
cx = 0; cy = 0;
% Compute charging field
Ch = arrayfun(@(xx, yy) (1/(1+exp(-t1*(xx-(cx-h))))).*(1/(1+exp(t1*(xx-(cx+h))))).*...
                (1/(1+exp(-t1*(yy-(cy-h))))).*(1/(1+exp(t1*(yy+(h-cy))))), Xg, Yg);

% Robot initial positions
d = 3; M = 3;
pos = zeros(M,2);
for m=1:M
    pos(m, :) = [m, 0];
end

% Plot
figure; hold on;
% Process plume surface
mesh(Xg, Yg, Z_true, 'FaceAlpha', 0.7, 'EdgeAlpha', 0.4);
% Charging station as semi-transparent surface offset slightly
surf(Xg, Yg, zmin + 0.1 + Ch*(zmax-zmin)*0.2, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.2);

% Robots as spheres
for m=1:M
    plot3(pos(m,1), pos(m,2), zmax*1.1, 'o', 'MarkerSize', 10, ...
          'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    text(pos(m,1), pos(m,2), zmax*1.2, sprintf('R%d', m), 'FontSize', 12);
end

% Axis labels and view
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Plume Height');
axis([xmin xmax ymin ymax zmin zmax*1.5]);
view([45 30]);
legend({'Plume','Charging Field','Robot Init.'}, 'Location', 'northeast');
title('Initial Map: Plume, Robots & Charger');
grid on;
