clear 
close all
clc

x = linspace(-10,10,100);
y = linspace(-10,10,100);

% 2) Turn them into a grid  
[X,Y] = meshgrid(x,y);

% 3) Compute the sigmoid in each direction, then multiply elementwise
Sx = 1 ./ (1 + exp(-2*X));  
Sy = 1 ./ (1 + exp(-2*Y));  
F  = Sx .* Sy;       % note the .*


% 4) Plot it as a surface  
figure
surf(X,Y,F)
shading interp      % smooth colors
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('Two-Dimensional Sigmoid Surface')

% y = linspace(-10,10,100);
% [X,Y] = meshgrid(x,y);
% 
% % element-wise logistic factors:
% Lx  = 1./(1 + exp(-X));
% Lx_neg = 1./(1 + exp( X));    % = 1 – Lx
% Ly  = 1./(1 + exp(-Y));
% Ly_neg = 1./(1 + exp( Y));    % = 1 – Ly
% 
% F = Lx .* Lx_neg .* Ly .* Ly_neg;  % element-wise multiply
% 
% surf(X, Y, F)
% xlabel('x'), ylabel('y'), zlabel('f(x,y)')
% shading interp   % smooth shading
% colorbar
