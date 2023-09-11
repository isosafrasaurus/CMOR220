% Pierce Zhang, CMOR220, Fall 2023, Project 1: Grow a fern
% fern.m
% Draws two ferns using matrix transformations on 2x1 point matrices set
% by probability thresholds, as specified by the Project 1 specifications.
% Last modified: 11 September 2023

% Driver function to execute part one and two code
function fern
%% PART 1: PLOTTING A SIMPLE FERN
simplefern(0.3994);
% reduce cutoff
% simplefern(0.2);
% increase cutoff
% simplefern(0.7);

%% ANSWER TO PART 1 QUESTION:
% When the cutoff is increased, the points tend to be more concentrated
% to the bottom of the fern graph and its branches. When the cutoff is
% decreased, the points are more concentrated towards the top.

%% PART 2: PLOTTING AN ADVANCED FERN
advancedfern();
end

% Inputs: float cutoff, a threshold for determining which form of the
% matrix transformation on z should be used through probability.
% Outputs: none
function simplefern(cutoff)
% This function will display a graph of 4000 points that resembles a
% fern by transforming the starting matrix z = [1 ; 0] based on two
% given transformations according to a probability set as cutoff for
% the first and (1 - cutoff) for the second.
figure()
hold on

% Title
title("Part 1: Simple fern, cutoff = " + cutoff)

% Starting point
z = [1 ; 0];

% Repeat for 4,000 iterations (doc specified)
for i = 1:4000
    % Plot point (z(1),z(2)) pointwise
    plot(z(1), z(2), '.');

    % Random number in [0,1]
    r = rand();
    if r < cutoff % Transformation 1
        z = [0.4 -0.3733 ; 0.06 0.6] * z + [0.3533 ; 0];
    else % Transformation 2
        z = [-0.8 -0.1867 ; 0.1371 0.80] * z + [1.1 ; 0.1];
    end
end
hold off
end

% Inputs: none
% Outputs: none
function advancedfern()
% This function will display a graph resembling a tilted fern by
% plotting all points z that are generated from a series of matrix
% transformations on z based on four probability thresholds.
figure()
hold on

% Title
title("Part 2: Advanced fern")

% Starting point
z = [1 ; 0];

% Repeat for 4,000 iterations (doc specified)
for i = 1:4000
    % Plot point (z(1),z(2)) pointwise
    plot(z(1), z(2), '.');

    % Random number in [0,1]
    r = rand();
    if r <= 0.01 % Transformation 1 (probability 0.01)
        z = [0 0 ; 0 0.16] * z;
    elseif r <= 0.76 % Transformation 2 (probability 0.75)
        z = [0.85 0.04 ; -0.04 0.85] * z + [0 ; 1.6];
    elseif r <= 0.88 % Transformation 3 (probability 0.12)
        z = [0.2 -0.26 ; 0.23 0.22] * z + [0 ; 1.6];
    else % Transformation 2 (probability 0.12)
        z = [-0.15 0.28 ; 0.26 0.24] * z + [0 ; 0.44];
    end
end
hold off
end