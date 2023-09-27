% Pierce Zhang, CMOR220, Fall 2023, Project 3B: Newton Fractals Bonus
% newton_bonus.m
% Bonus for Newton fractals project
% Last modified: 27 September 2023

% Driver function
% Uses qnewt to produce four plots for each of four quartic functions
function newton_bonus
% Parameters of Newton's method
xt = [-0.2,0.0025,0.2];
yt = [-0.2,0.0025,0.2];
maxiter = 20;

qnewt([2 -4 2-1i -1i 100],xt,yt,maxiter);
end

% Inputs: q, a vector of coefficients of a polynomial
% Outputs: dq, a vector of coefficients of the analytical derivative of q
function [dq] = myownpolyder(q)
% Evaluates the analytical derivative of a polynomial
% Derivative function should have one less term
dq = zeros(1,length(q) - 1);
for i = 1:length(dq)
    % Where the coefficient is q(i) and the power is length(q) - i
    dq(i) = q(i) * (length(q) - i);
end
end

% Inputs:
% coeff, the coefficient of a monomial term
% img, a boolean indicating whether the coefficient is imaginary
% power, the power of z in the monomial term
% Outputs: term, a string containing the formatted monomial term
function [term] = qterm(coeff, img, power)
term = "";
% Sign handler
if (coeff < 0)
    term = term + '-';
elseif (coeff > 0)
    term = term + '+';
else
    return
end
% Coefficient handler
if (abs(coeff) ~= 1)
    term = term + abs(coeff);
elseif (power == 0 && img == false)
    term = term + abs(coeff);
end
%Imaginary unit handler
if (img)
    term = term + 'i';
end
% Power handler
if (power == 1)
    term = term + 'z';
elseif (power > 1)
    term = term + 'z^' + power;
end
end

% Inputs: q, a  vector of coefficients for a polynomial
% Outputs: name, a string which represents a polynomial
function [name] = qlab(q)
% Produces a title for a given polynomial as a coefficient vector
name = "";
for i=1:length(q)
    coeff = q(i);
    power = length(q) - i;
    % Skip over all 0 coefficient terms
    if (coeff ~= 0)
        name = name + qterm(real(coeff), false, power);
        name = name + qterm(imag(coeff), true, power);
    end
end
% Remove leading plus sign.
if (name ~= "" && extract(name,1) == '+')
    name = eraseBetween(name, 1, 1);
end
end

% Inputs: 
% Outputs: color, verdict based on
function [color] = newt(q, dq, roots, z, maxiter)
rootcolors = ['r','g','b','y'];
speccolors = ['k','m',"#D95319"];

color = "";
for i=1:maxiter
    qval = polyval(q,z);
    dqval = polyval(dq, z);

    if (dqval == 0)
        color = speccolors(1);
        break
    end
    if (z > 1e9)
        color = speccolors(2);
        break
    end
    z = z - qval/dqval;

    for ri=1:length(roots)
        if (abs(z - roots(ri))<0.001)
            color = rootcolors(ri);
            break
        end
    end
end
end

% Inputs:
% q, a vector of up to 5 complex coefficients of a polynomial
% xt, a vector of 3 x grid values, xlo, xinc, and xhi
% yt, a vector of 3 y grid values, ylo, yinc, and yhi
% maxiter, the maximum number of iterations
% Outputs: none (produces plots)
function qnewt(q, xt, yt, maxiter)
% Plots the Newton's basins and wasteland for a given function
% 1. Generate initial guesses using meshgrid,xt,yt
% 2. Apply Newton's method
% 3. Assign the colors by using find command
% 4. Apply qlab to assign the titles for the plots

for x = xt(1):xt(2):xt(3)
    for y = yt(1):yt(2):yt(3)
        color = newt(q, myownpolyder(q), roots(q), x+1i*y,maxiter);
        if (color ~= "")
            plot(x,y,'.','Color',color,'MarkerSize',5)
        end
        hold on
    end
end
title(qlab(q));
end