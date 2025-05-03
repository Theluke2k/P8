clear
close all
clc

import casadi.*

% Create an opti stack
opti = Opti();

% Create decision variables
x = opti.variable(2,1);

% Create parameters
p = opti.parameter(2,1);

% Define constraints
opti.subject_to(x >= p)

% Define objective function
opti.minimize( x(1,1) + x(2,1) );

% Configure printing for ipopt solver
opts = struct();
opts.ipopt.print_level = false;

opti.set_value(p,[5; 11]);

% Choose the ipopt solver
opti.solver('ipopt', opts);
sol = opti.solve();
x_opt= sol.value(x);
disp(x_opt)
