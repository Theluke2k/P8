clear
close all
clc

import casadi.*

% Create an opti stack
opti = Opti();

% Create decision variables
%x1 = opti.variable();
x2 = opti.variable();

% Create parameters
p = opti.parameter();

% Define constraints
%opti.subject_to(x2 - x1 == p)
%x1 = x2-p;
opti.subject_to(x2-p+x2 >= 0)

% Define objective function
opti.minimize( x2-p +  x2);

% Configure printing for ipopt solver
opts = struct();
opts.ipopt.print_level = false;

opti.set_value(p,5);

% Choose the ipopt solver
opti.solver('ipopt', opts);
sol = opti.solve();
%x1_opt= sol.value(x1);
x2_opt= sol.value(x2);
%disp(x1_opt)
disp(x2_opt)
% % Solve for different values of p
% for p_val = [2, -1]
%     opti.set_value(p, p_val);
%     sol = opti.solve();
%     x_opt= sol.value(x);
%     fprintf('For p = %g  →  x* = %g\n', p_val, x_opt);
% end