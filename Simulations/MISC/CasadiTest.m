clear
close all
clc

import casadi.*

% Create an opti stack
opti = Opti();

% Create decision variables
x = opti.variable();

% Create parameters
%p = opti.parameter();

% Define objective function
opti.minimize( (x-2)^2 );

% Define constraints
opti.subject_to( x >= 0);

% Configure printing for ipopt solver
opts = struct();
opts.ipopt.print_level = false;

% Choose the ipopt solver
opti.solver('ipopt', opts);
sol = opti.solve();
x_opt= sol.value(x);
disp(x_opt)
% % Solve for different values of p
% for p_val = [2, -1]
%     opti.set_value(p, p_val);
%     sol = opti.solve();
%     x_opt= sol.value(x);
%     fprintf('For p = %g  →  x* = %g\n', p_val, x_opt);
% end