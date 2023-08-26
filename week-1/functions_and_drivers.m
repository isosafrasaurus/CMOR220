% Pierce Zhang, CMOR220, FALL 2023, Competency functions and drivers
% functions_and_drivers.m
% Answers to functions and drivers competency questions
% Last modified: August 26, 2023

% Driver to answer the questions
function functions_and_drivers
    % Problem 1
    [y1] = harm(1,2)
    [y2] = harm(2,4)
    [y3] = harm(3,6)

    % Problem 2
    [sinout1, cosout1, tanout1] = trigs(pi / 3)
    [sinout2, cosout2, tanout2] = trigs(pi / 6)
    
    % Problem 3
    [x1, x2] = quad(1/3, 5, 3)

    % Problem 4
    H = 0.00344; A = 0.00344; HA = 0.66;
    [Ka, pKa, G] = equil(H, A, HA)

    % Problem 5
    lineqs(-2, 8, 3, 6);
end

% Inputs: a, the amplitude, and t, the time of the simple harmonic motion
% equation
% Outputs: y, the value of the function given a and t
function [y] = harm(a, t)
    % This function will output the value of a simple harmonic motion
    % function for a given amplitude and time, with angular frequency set
    % to 3 and y-offset to +4.
    y = a * sin(3 * t) + 4;
end

% Inputs: x, the input angle in radians
% Outputs: sinout, the value of sin(x), cosout, the value of cos(x), and
% tanout, the value of tan(x)
function [sinout, cosout, tanout] = trigs(x)
    sinout = sin(x);
    cosout = cos(x);
    tanout = tan(x);
end

% Inputs: coefficients a, b, and c, of a second-degree polynomial
% Outputs: x1 and x2, the roots of said polynomial
function [x1, x2] = quad(a, b, c)
    det = sqrt(b^2 - 4*a*c);
    x1 = (-b + det) / (2 * a);
    x2 = (-b - det) / (2 * a);
end

% Inputs: H, the molar concentration of protons in solution, A, the molar
% concentration of CH3OO-, and HA, the molar concentration of acetic acid
% (all prior to reaction).
% Outputs: Ka, the value of the equilibrium constant, pKa, the "power" of
% the equilibrium constant, and G, the Gibbs free energy of dissociation
function [Ka, pKa, G] = equil(H, A, HA)
    R = 8.1345; % units J / K
    T = 298;

    Ka = H * A / HA;
    pKa = -log10(Ka);
    G = -R * T * log(Ka);
end

% Inputs: x1, y1, the coordinates of the first point, and x2, y2, the
% coordinates of the second
% Outputs: (function does not return a value) a printout of the equations
% of the line that passes through (x1, y1) and (x2, y2) in both point-slope
% and slope-intercept form.
function lineqs(x1, y1, x2, y2)
    m = (y2 - y1) / (x2 - x1);
    disp("y - (" + y1 + ") = (" + m + ")(x - (" + x1 + "))");
    disp("y = (" + m + ")x + (" + (-m * x1 + y1) + ")");
end