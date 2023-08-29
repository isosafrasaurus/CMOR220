% Pierce Zhang, CMOR220, FALL 2023, Compentency if-statements and code loops
% if_statements_and_loops.m
% Answers to if-statements and code loops competency
% Last modified: August 28, 2023

% Driver to answer the questions
function if_statements_and_loops
    % Problem 1
    display_power(9);
    display_power(99);
    display_power(999);
    display_power(1000);

    % Problem 2
    Problem2Result = sinexpthing(10)

    % Problem 3
    Problem3Result_1 = summation(6.7, 125)
    Problem3Result_2 = summation(2.55, 298)
    Problem3Result_3 = summation(5.07, 373)

    % Problem 4
    gamble(50, 50);
    gamble(50, 150);
    gamble(50, 1000);

    % Problem 5
    HelloGoodbye1or2(1);
    HelloGoodbyeN1to2(1);
    HelloGoodbyene1(1);
end

% Inputs: n, the number to be checked
% Outputs: [will display to console] "Value is less than {power ceil of n}"
function display_power(n)
    % This function will conditionally check and display "Value is less 
    % than {power ceil of n}" for given inputs of n, going up to 10^3.
    if (n < 10^1)
        disp("Value is less than 10^1");
    elseif (n < 10^2)
        disp("Value is less than 10^2");
    elseif (n < 10^3)
        disp("Value is less than 10^3");
    else
        disp("Woah");
    end
end

% Inputs: n, the number of times the function should be recursed
% Outputs: f, the value of f(n) = sin(exp(1.01 - f(n-1)))
function [fn] = sinexpthing(n)
    % Use for-loop to recurse
    f = sin(exp(1.01)); % f(0)
    for i=1:n
        f = sin(exp(1.01 - f));
    end
    fn = f;
end

% Inputs: x, the value to be added incrementally, and T, the threshold
% Outputs: t, the number of times X must be added for sum to exceed T
function [t] = summation(x, T)
    t = 0;
    sum = 0;
    while(sum <= T)
        sum = sum + x;
        t = t + 1;
    end
end

% Inputs: x0, the principal, and N, the number of times to play
% Outputs: [will print to console] x, the new total
function gamble(x0, N)
    x = x0;
    for n=1:N
        r1 = rand(); r2 = rand();
        if (n <= 100)
            x=x+r1*x0-r2*x0;
        elseif (n <= 500)
            x=x+(r1-0.1)*x0-(r2+0.125)*x0;
        else
            x=x+(r1-0.2)*x0-(r2+0.25)*x0;
        end
    end
    disp(x);
end

% Inputs: input
% Outputs: [will print to console] "Hello!" whenever the input is 1 or 2, and Goodbye! otherwise
function HelloGoodbye1or2(input)
    if (input == 1 || input == 2)
        disp("Hello!");
    else 
        disp("Goodbye!");
    end
end

% Inputs: input
% Outputs: "Hello!" for any input strictly greater than -1 and strictly less than 2, and "Goodbye!" otherwise
function HelloGoodbyeN1to2(input)
    if (input > -1 && input < 2)
        disp("Hello!");
    else
        disp("Goodbye!");
    end
end

% Inputs: input
% Outputs: "Hello!" for any input not equal to 1, and "Goodbye!" otherwise
function HelloGoodbyene1(input)
    if (input ~= 1)
        disp("Hello!")
    else
        disp("Goodbye!")
    end
end