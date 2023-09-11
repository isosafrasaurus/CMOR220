% Pierce Zhang, CMOR220, FALL 2023, Compentency on fern project
% competency_fern_project.m
% Answers to fern project competency
% Last modified: 9 September 2023

% Driver to answer the questions
function competency_fern_project
    % Problem 1
    game1 = problem1(50)
    game2 = problem1(100)
    game3 = problem1(500)
    % Answer: NO. You shouldn't play.

    % Problem 2
    problem2()

    % Problem 3
    problem3()
end

% Inputs: x0, the principal
% Outputs: x, the ending amount of money after 500 games
function [x] = problem1(x0)
    % Description: performs a for-loop with interal rand in order to
    % estimate the gambling outcome as indicated by the specs and
    % input/output.
    x = x0;
    for i = 1:500
        r = rand();
        if r <= 0.3
            x = x + 1.15 * x0;
        elseif r <= 0.80
            x = x - 0.50 * x0;
        else
            x = x - 0.90 * x0;
        end
    end
end

% Inputs: none
% Outputs: none
function problem2()
    % Description: displays the plot of y = sqrt(x) for x in [0, 10] using
    % 200 points. There is a 10% chance the points displayed are red, 30%
    % they are blue, and 60% they are green.
    figure()
    hold on
    x = linspace(0,10,200);
    for i = x
        r = rand();
        if (r <= 0.1)
            plot(i, sqrt(i), ".",'Color', 'r')
        elseif (r <= 0.40)
            plot(i, sqrt(i), ".",'Color', 'b')
        else
            plot(i, sqrt(i), ".",'Color', 'g')
        end
    end
    hold off
end

% Inputs: none
% Outputs: none
function problem3()
    % Description: Displays a fractal over 1000 point iterations using a
    % specified vector transformation on the starting vector z = [1;0].
    figure()
    hold on
    z = [1;0];
    theta = pi / 2.5;
    for i = 1:1000
        plot(z(1), z(2), '.')
        z = [cos(theta)+0.025, sin(theta); -sin(theta), cos(theta)]*z;  
    end
    hold off
end