function [h] = get_h(states,px,py,sc)
    h = (sc(1)*states(1)*sc(3)*states(3)/pi)*exp(-sc(3)*states(3)*((px-sc(5)*states(5)).^2 + (py-sc(7)*states(7)).^2));
end

