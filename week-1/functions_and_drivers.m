% Pierce Zhang, CMOR220, Fall 2023, Competency functions and drivers
% functions_and_drivers.m
% Answers to functions and drivers competency
% 24 August 2023

function functions_and_drivers
    % Problem 1
    disp(harm(1,2));
    disp(harm(2,4));
    disp(harm(3,6));

    % Problem 2
    [sinout, cosout, tanout] = trigs(pi / 3);
    disp(sinout); disp(cosout); disp(tanout);
    
    % Problem 3
    a = 1/3; b = 5; c = 3;
    [x1, x2] = quad(a, b, c);

    disp(x1);
    disp(x2);

    % Problem 4
    disp("EQUILIBRIUM OUTPUT WILL FOLLOW");

    H = 0.00344; A = 0.00344; HA = 0.66;
    [Ka, pKa, G] = equil(H, A, HA);

    disp(Ka + " " + pKa + " " + G);

    % Problem 5
    lineqs(-2, 8, 3, 6);
end

% Inputs:
% Outputs:
function [y] = harm(a, t)
    y = a * sin(3 * t) + 4;
end

function [sinout, cosout, tanout] = trigs(x)
    sinout = sin(x);
    cosout = cos(x);
    tanout = tan(x);
end

function [x1, x2] = quad(a, b, c)
    det = sqrt(b^2 - 4*a*c);
    x1 = (-b + det) / (2 * a);
    x2 = (-b - det) / (2 * a);
end

function [Ka, pKa, G] = equil(H, A, HA)
    R = 8.1345; % J / K
    T = 298;

    Ka = H * A / HA;
    pKa = -log10(Ka);
    G = -R * T * log(Ka);
end

function lineqs(x1, y1, x2, y2)
    m = (y2 - y1) / (x2 - x1);
    disp("y - (" + y1 + ") = (" + m + ")(x - (" + x1 + "))");
    disp("y = (" + m + ")x + (" + (-m * x1 + y1) + ")");
end