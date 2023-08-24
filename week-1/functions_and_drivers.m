function functions_and_drivers
    disp("p. 1 where a = 1, t = 2: " + harm(1,2));
    disp("p. 1 where a = 2, t = 4: " + harm(2,4));
    disp("p. 1 where a = 3, t = 6: " + harm(3,6));

    
end

function [y] = harm(a, t)
    y = a * sin(3 * t) + 4;
end

function [sinout, cosout, tanout] = trigs(x)
    sinout = sin(x);
    cosout = cos(x);
    tanout = tan(x);
end