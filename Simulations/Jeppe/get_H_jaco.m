function [H] = get_H_jaco(z, x, M,sc)
    Nx_r = length(z);   % Number of states
    for m = 1:M
        px = x(2*m-1);  % Robot x-position
        py = x(2*m);    % Robot y_position

        % Build row of jacobian
        row = zeros(1,Nx_r);
        row(1,1) = get_h(z, px, py, sc)/z(1);
        row(1,3) = get_h(z, px, py, sc)*(1/z(3) - sc(3)*((px - sc(5)*z(5))^2 + (py - sc(7)*z(7))^2));
        row(1,5) = get_h(z, px, py, sc)*(-2*sc(5)*sc(3)*z(3)*(sc(5)*z(5) - px));
        row(1,7) = get_h(z, px, py, sc)*(-2*sc(7)*sc(3)*z(3)*(sc(7)*z(7) - py));

        % Insert row
        H(m, :) = row;
    end
end