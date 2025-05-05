clear
close all
clc

import casadi.*

% Create an opti stack
opti = Opti();

% Create decision variables
du = opti.variable();
x = opti.variable();

du = 

% Define constraints
opti.subject_to( x+z >= 0)

% Define objective function
opti.minimize( y+z );

% Configure printing for ipopt solver
opts = struct();
opts.ipopt.print_level = false;

% Choose the ipopt solver
opti.solver('ipopt', opts);
sol = opti.solve();
sol.value(y)
sol.value(z)