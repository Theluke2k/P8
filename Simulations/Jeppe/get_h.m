function [h] = get_h(states,px,py)
    h = (states(1)*states(3)/pi)*exp(-states(3)*((px-states(5)).^2 + (py-states(7)).^2));
end

