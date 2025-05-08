clear
close all
clc

import casadi.*

% Create an opti stack
opti = Opti();

% Create decision variables
x = opti.variable(3,3);
y = MX.zeros(3,3);
z = 0;
for i = 1:3
    y(:,i) = i*x(:,i);
    z = z + y(i,i);
end

% y = [1*x11 2*x12 2*x13
%      1*x21 2*x22 2*x23
%      1*x31 2*x32 2*x33]

opti.subject_to(x(:,1) >= 1);
opti.subject_to(x(:,2) >= 2);
opti.subject_to(x(:,3) >= 3);

% Define objective function
opti.minimize( z );

% Configure printing for ipopt solver
opts = struct();
opts.ipopt.print_level = false;

% Choose the ipopt solver
opti.solver('ipopt', opts);

sol = opti.solve()

sol.value(z)
sol.value(x)