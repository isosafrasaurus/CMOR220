% Pierce Zhang, CMOR220, Fall 2023, Competency functions and drivers
% functions_and_drivers.m
% Answers to functions and drivers competency
% 24 August 2023

% Driver
function functions_and_drivers
    % Answer to problem 1
    disp("p. 1 where a = 1, t = 2: " + Question1(1,2));
    disp("p. 1 where a = 2, t = 4: " + Question1(2,4));
    disp("p. 1 where a = 3, t = 6: " + Question1(3,6));    
end

% Inputs:
% Outputs:
function [y] = Question1(a, t)
    y = a * sin(3 * t) + 4;
end

function [sinout, cosout, tanout] = Question2(x)
    sinout = sin(x);
    cosout = cos(x);
    tanout = tan(x);
end