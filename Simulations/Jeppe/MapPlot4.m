%% ——————————————————————————————————————————————
%  Figure defaults: LaTeX interpreter & CMU Serif
%  Place this at the top of your script, before any figure()
% ——————————————————————————————————————————————
% Figure defaults (landscape 19×12 cm)
set(groot, ...
  'defaultFigureUnits','centimeters', ...
  'defaultFigurePosition',[2 2 25 11], ...       % 19 cm wide, 12 cm tall
  'defaultAxesFontName','CMU Serif', ...
  'defaultTextFontName','CMU Serif', ...
  'defaultAxesTickLabelInterpreter','latex', ...
  'defaultTextInterpreter','latex', ...
  'defaultLegendInterpreter','latex', ...
  'defaultAxesFontSize',10, ...
  'defaultTextFontSize',10, ...
  'defaultLegendFontSize',9);


%% ——————————————————————————————————————————————
clear; close all; clc;

% Map bounds and grid
xmin = 0; xmax = 40;
ymin = 0; ymax = 40;
zmin = 0; zmax = 5;
x_axis = linspace(xmin, xmax, 100);
y_axis = linspace(ymin, ymax, 100)';
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
    (1/(1+exp( t1*(xx-(cx+h_ch))))) * ...
    (1/(1+exp(-t1*(yy-(cy-h_ch))))) * ...
    (1/(1+exp( t1*(yy-(cy+h_ch))))) * charge_rate, Xg, Yg);

% Robot initials and colors
M = 3;
pos = [(1:M)', zeros(M,1)];
colors = [0.1216 0.4667 0.7059; 1.0000 0.4980 0.0549; 0.5804 0.4039 0.7412];
thickness = 7;  % original marker size

% Create figure
hFig = figure('Color','w');
lightangle(45,30); lightangle(-30,60);

% --- Left: Plume + Robots ---
ax1 = subplot(1,2,1); hold(ax1,'on');
surf(ax1, Xg, Yg, Z_true, ...
    'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.1,'FaceAlpha',0.9);
contour3(ax1, Xg, Yg, Z_true, linspace(0,max(Z_true(:)),8),...
    'k','LineWidth',0.5);
for m=1:M
    plot3(ax1, pos(m,1), pos(m,2), 0, 'o', ...
        'MarkerSize', thickness, ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end
colormap(ax1, parula);shading(ax1,'interp');
lighting(ax1,'gouraud'); material(ax1,'dull');
axis(ax1,[xmin xmax ymin ymax 0 4]);
view(ax1,45,30);
xlabel(ax1,'$x$~[m]');
ylabel(ax1,'$y$~[m]');
zlabel(ax1,'Process~Amplitude');
title(ax1,'\textbf{Initial Process and Robots}');
grid(ax1,'on'); box(ax1,'on');

% --- Right: Charger + Robots ---
ax2 = subplot(1,2,2); hold(ax2,'on');
surf(ax2, Xg, Yg, Ch*5, ...
    'FaceColor','interp','EdgeColor','k','EdgeAlpha',0.1,'FaceAlpha',0.7);
contour3(ax2, Xg, Yg, Ch*5, ...
    linspace(min(Ch(:))*5,max(Ch(:))*5,6),'k','LineWidth',0.5);
for m=1:M
    plot3(ax2, pos(m,1), pos(m,2), 0, 'o', ...
        'MarkerSize', thickness, ...
        'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end
colormap(ax2, autumn);shading(ax2,'interp');
lighting(ax2,'gouraud'); material(ax2,'dull');
axis(ax2,[xmin xmax ymin ymax 0 2]);
view(ax2,45,30);
xlabel(ax2,'$x$ [m]');
ylabel(ax2,'$y$ [m]');
zlabel(ax2,'Charge~Rate');
title(ax2,'\textbf{Charging Station and Robots}');
grid(ax2,'on'); box(ax2,'on');

% size the figure window first
set(hFig,'Units','centimeters','Position',[2 2 25 11]);

% then export
exportgraphics(hFig, 'ProcessAndChargingPlot.png', ...
               'Resolution', 2000);   % dpi
