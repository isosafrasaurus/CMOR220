% Pierce Zhang, CMOR220, FALL 2023, Competency Euler's Method and ODEs
% competency_eulers_method.m
% Answers to competency on Euler' method and ODEs
% Last modified: 4 October 2023

% Driver
function competency_eulers_method
    [y] = eulers_method(@prob1_fun,0,4,40,0.2);
    figure; hold on; grid on;
    plot(0:0.2:40,y);

    clear;

    [y] = eulers_method(@prob1_fun,0,4,40,0.01);
    figure; hold on; grid on;
    plot(0:0.01:40,y);

    clear;

    [x,y] = eulers_method_sys(@prob3_funx,@prob3_funy,1,0,0,40,0.01);
    figure; hold on; grid on;
    plot(x,y);

    clear;
        
    dTdt = @(t, T) -k*(T(t) - TA);
end

function [yprime] = prob1_fun(x, y)
    yprime = 0.05 * x * y * sin(x) * cos(2.5*x);
end

function [y] = eulers_method(fun,x0,y0,b,dx)
    x = x0:dx:b;
    y=zeros(1,length(x));
    y(1)=y0;
    for n=1:length(x)-1
        dy = (fun(x(n),y(n)))*dx;
        y(n+1)=y(n)+dy;
    end
end

function [xprime] = prob3_funx(x,y)
    xprime = y + 0.2*x;
end

function [yprime] = prob3_funy(x,y)
    yprime = (1 - x^2)*y - x;
end

function [x,y] = eulers_method_sys(funx,funy,x0,y0,t,tf,dt)
    tdomain = t:dt:tf;
    x=zeros(1,length(tdomain));
    y=zeros(1,length(tdomain));
    x(1)=x0;
    y(1)=y0;
    for n=1:length(tdomain)-1
        dx = (funx(x(n),y(n)))*dt;
        dy = (funy(x(n),y(n)))*dt;
        x(n+1)=x(n)+dx;
        y(n+1)=y(n)+dy;
    end
end