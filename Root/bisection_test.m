function bisection_test
    root = bisection_method(1)
end

function root = bisection_method(L)
    % Define the function
    f = @(x) sin(x*L) + x*cos(x*L);
    
    % Define the initial interval [a, b]
    a = 0.1;
    b = 3;
    
    % Define the tolerance
    tol = 1e-6;
    
    % Perform the bisection method
    while abs(b - a) > tol
        c = (a + b) / 2;
        if f(a) * f(c) < 0
            b = c;
        else
            a = c;
        end
    end
    
    % The root is the midpoint of the final interval
    root = (a + b) / 2;
end