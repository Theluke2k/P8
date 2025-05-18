clear
close all
clc

%% Initial Conditions Combined Map Plot
% Overlaying plume, charging field, and robots in one 3D visualization
clear; close all; clc;

% Map bounds and grid
dim = struct('xmin',0,'xmax',40,'ymin',0,'ymax',40,'zmin',0,'zmax',5);
x_axis = linspace(dim.xmin, dim.xmax, 150);
y_axis = linspace(dim.ymin, dim.ymax, 150)';
[Xg, Yg] = meshgrid(x_axis, y_axis);

% Process initial state and scaling
sc = [200 1 0.05 0.0001 20 0.1 20 0.1]';
T = diag(sc);
z0 = [200;0;0.05;0;20;0;20;0];
z0s = inv(T)*z0;
Z_true = get_h(z0s, x_axis, y_axis, sc);

% Charging station field
t1 = 0.5; h_ch = 5; cx = 0; cy = 0; charge_rate = 0.2;
Ch = arrayfun(@(xx,yy) ...
    (1/(1+exp(-t1*(xx-(cx-h_ch))))) * ...
    (1/(1+exp(t1*(xx+(h_ch-cx))))) * ...
    (1/(1+exp(-t1*(yy-(cy-h_ch))))) * ...
    (1/(1+exp(t1*(yy+(h_ch-cy))))) * charge_rate, Xg, Yg);

% Robots initial positions and colors
M = 3;
pos = [(1:M)', zeros(M,1)];
thickness = 5;
colors = [ ...
   0.1216   0.4667   0.7059;  % blue
   1.0000   0.4980   0.0549;  % orange
   0.5804   0.4039   0.7412]; % purple
robot_z = ones(M,1)*1;  % z = 1 for all robots

% Create figure
figure('Color','w','Position',[200 200 800 600]); hold on;

% --- Plume surface ---
hP = surf(Xg, Yg, Z_true, ...
    'FaceColor','interp','EdgeColor','none','FaceAlpha',0.8);
colormap(parula);
contour3(Xg, Yg, Z_true, 8, 'k', 'LineWidth',0.4);

% --- Charging field surface at ground level ---
hC = surf(Xg, Yg, Ch, ...
    'FaceColor','interp','EdgeColor','none','FaceAlpha',0.6);
contour3(Xg, Yg, Ch, 6, 'k', 'LineWidth',0.4);

% --- Robots (matching original style) ---
for m = 1:M
    plot3(pos(m,1), pos(m,2), robot_z(m), 'o', ...
        'MarkerSize', thickness, ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end

% --- Styling ---
shading interp; lighting gouraud; material dull;
axis([dim.xmin dim.xmax dim.ymin dim.ymax dim.zmin dim.zmax]);
view(45,30);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Value');
title('Initial Plume + Charging Field + Robots');
grid on; box on;

% Colorbars for each field
tbb = colorbar; tbb.Label.String = 'Plume Height';
cb2 = colorbar('Location','southoutside'); cb2.Label.String = 'Charge Rate';

hold off;
