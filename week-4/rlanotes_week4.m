function rlanotes_week4
    f = @(x) x^4 + x^2 + 16
    a = 0;
    b = pi;
    tol = 1e-6;
    [root, iterations] = bisection_method(f, a, b, tol);
    disp(root + " " + iterations)

function [root, iterations] = bisection_method(f, a, b, tol)
    
    iterations = ceil(log2((b - a) / tol));
    for i = 1:iterations
        c = (a + b) / 2;
        if f(c) == 0 || (b - a) / 2 < tol
            root = c;
            return;     
        end
        if sign(f(c)) == sign(f(a))
            a = c;
        else
            b = c;
        end
    end
    
    root = (a + b) / 2;