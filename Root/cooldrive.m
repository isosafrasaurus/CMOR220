function cooldrive()
    % set a, t and a range of L and b and find x
    % and plot x^2 against L
    a = 0.1; b = 3; t = 0.01;
    L = 1 : 0.1 : 4;

    x = zeros(0,length(L));
    for cnt = 1:length(L)
        [x(cnt), iter] = debis(a,b,t,L(cnt));
        b = b - 0.05;
    end
    
    figure()
    plot(L, x.^2, "-x")
    title("Cooling Rate vs. Bar Length")
    xlabel("Length, L")
    ylabel("Decay Rate, x^2")

    % End of Problem 1

    clear

    % Problem 2

    L = 1; a = 0.1; b = 3; x0 = (a+b)/2;
    t = 10.^-(1:8);

    x_bis = zeros(0,length(t));
    x_newt = zeros(0,length(t));
    for cnt = 1:length(t)
        [~,x_bis(cnt)] = debis(a,b,t(cnt),L);
        [~,x_newt(cnt)] = denewt(x0,t(cnt),L);
    end

    figure()
    semilogx(t, x_newt, "-o","Color",'r');
    hold on
    semilogx(t, x_bis, "-x","Color",'b');
    hold on
    title("Bisection vs. Newton")
    xlabel("Tolerance")
    ylabel("# of iterations")
    legend("Newton","Bisection")

end


function[x,iter] = denewt(x,t,L)
    %Newtons method, calls coolfun and coolfundx
    iter = 0;
    while (abs(coolfun(x,L)) > t)
        x = x - coolfun(x,L)/coolfundx(x,L);
        iter = iter + 1;
    end
end

function [x,iter] = debis(a,b,tol,L)
    % Perform the bisection method
    iter = 0;
    while abs(b - a) > tol
        x = (a + b) / 2;
        if coolfun(a,L) * coolfun(x,L) < 0
            b = x;
        else
            a = x;
        end
        iter = iter + 1;
    end
    
    x = (a + b) / 2;
end

function val = coolfun(x,L) %for a given x and L evaluate "cool"
    val = sin(x*L) + x*cos(x*L);
end

function val = coolfundx(x,L)
    % evaluate the derivative, with respect to x, of coolfun
    val = L*cos(x*L) + cos(x*L) - x*L*sin(x*L);
end