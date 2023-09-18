%Inthesamedriver
%seta,tandarangeofLandbandfindx
%andplotx^2againstL

%callsdebisanddenewtandplotstheiriterspertol

function cooldrive()
    % set a, t and a range of L and b and find x
    % and plot x^2 against L
    a = 0.1; b = 0.3; t = 0.01;
    L = 1 : 0.1 : 4;

    x = zeros(length(L));
    for cnt = 1:length(L)
        x(cnt) = debis(a,b,t,L(cnt));
    end
    plot(L, x.^2)
end


function[x,iter]=denewt(x,tol,L)
    %Newtons method, calls coolfun and coolfundx
    while (abs(coolfun(x,L)) > tol)
        x = x - coolfun(x,L)/coolfundx(x,L);
    end
end

function[x,iter]=debis(a,b,t,L) % for a given a, b, t, and L find x
    % code the bisection method here
    x = a;
    while (abs(coolfun(x,L)) > t)
        x = (a + b) / 2;
        if (coolfun(a,L) * coolfun(b,L) > 0)
            disp("no solution in a,b");
        elseif coolfun(a,L) * coolfun(x,L) < 0
            b = x;
        else
            a = x;
        end
    end
end

function val = coolfun(x,L) %for a given x and L evaluate "cool"
    val = sin(x*L) + x*cos(x*L);
end

function val = coolfundx(x,L)
    % evaluate the derivative, with respect to x, of coolfun
    val = L*cos(x*L) + cos(x*L) - x*L*sin(x*L);
end