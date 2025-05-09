function [h_vec] = get_h_vec(z,x,M,sc)
    h_vec = [];
    for m = 1:M
        h_vec = [h_vec; get_h(z, x(2*m-1), x(2*m),sc)];
    end
end

