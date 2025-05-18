clear
close all
clc

%% Initial Conditions Map Plot
% Two styled subplots: plume+robots and charger+robots with subtle texture
clear; close all; clc;

% Map bounds and grid
xmin = 0; xmax = 40;
ymin = 0; ymax = 40;
zmin = 0; zmax = 5;
[x_axis, y_axis] = deal(linspace(xmin,xmax,100), linspace(ymin,ymax,100)');
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
    (1/(1+exp(t1*(xx-(cx+h_ch))))) * ...
    (1/(1+exp(-t1*(yy-(cy-h_ch))))) * ...
    (1/(1+exp(t1*(yy-(cy+h_ch))))) * charge_rate, Xg, Yg);

% Robot initials
M = 3; pos = [(1:M)', zeros(M,1)];

% Setup figure
gcf = figure('Position',[100 100 1200 500], 'Color','w');

% Lighting for both plots
lightangle(45,30); lightangle(-30,60);

% --- Left: Plume + Robots ---
ax1 = subplot(1,2,1); hold(ax1,'on');
hp = surf(ax1, Xg, Yg, Z_true, ...
    'FaceColor','interp', 'EdgeColor','k', 'EdgeAlpha',0.1, 'FaceAlpha',0.9);
% add contour for texture
contourLevels = linspace(0,max(Z_true(:)),8);
[c,hc] = contour3(ax1, Xg, Yg, Z_true, contourLevels, 'k');
set(hc,'LineWidth',0.5,'LineStyle','-');

% robots
plot3(ax1, pos(:,1), pos(:,2), 1*ones(M,1), 'o', ...
    'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','w','LineWidth',1.5);

% styling
colormap(ax1, parula);
shading(ax1,'interp');
lighting(ax1,'gouraud'); material(ax1,'dull');
axis(ax1,[xmin xmax ymin ymax zmin zmax]);
view(ax1,45,30);
xlabel(ax1,'X [m]'); ylabel(ax1,'Y [m]'); zlabel(ax1,'Plume Height');
title(ax1,'Initial Plume & Robots');
grid(ax1,'off'); box(ax1,'on');

% --- Right: Charger + Robots ---
ax2 = subplot(1,2,2); hold(ax2,'on');
hc2 = surf(ax2, Xg, Yg, Ch*5, ...
    'FaceColor','interp', 'EdgeColor','k', 'EdgeAlpha',0.1, 'FaceAlpha',0.7);
% subtle contour
contourLevels2 = linspace(min(Ch(:))*5, max(Ch(:))*5, 6);
[c2,hc2c] = contour3(ax2, Xg, Yg, Ch*5, contourLevels2, 'k');
set(hc2c,'LineWidth',0.5,'LineStyle','-');

plot3(ax2, pos(:,1), pos(:,2), 1*ones(M,1), 'o', ...
    'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','w','LineWidth',1.5);

colormap(ax2, autumn);
shading(ax2,'interp');
lighting(ax2,'gouraud'); material(ax2,'dull');
axis(ax2,[xmin xmax ymin ymax zmin zmax*1.2]);
view(ax2,45,30);
xlabel(ax2,'X [m]'); ylabel(ax2,'Y [m]'); zlabel(ax2,'Charge Rate');
title(ax2,'Initial Charging Field & Robots');
grid(ax2,'on'); box(ax2,'on');
grid(ax1,'on');