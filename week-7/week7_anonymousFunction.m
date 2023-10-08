% Anastasiya Protasov, CAAM 210, Fall 2023, anonymous functions
% week7_anonymousFunction.m
% Script to give examples of using function handles
% list of input: none
% list of output: none
% Last Modified: October 2, 2023

%% if you vectorize, your function can handle a variety of input


pts = linspace(-2,2,5); % define a vector

g = @(x) 3*x.^2;
g(pts)

[X,Y] = meshgrid(pts,pts);
ptMatrix = X+Y
g(ptMatrix)
 
%% function handles can make plotting curves easy

plot(pts,g(pts))
xlabel('x'); ylabel('y')
title('the parabola y = 3x^2')

%% function handles can have multiple variables

f2 = @(x,y) x.^2 + y.^2;
f2(pts,pts)
%  
surf(X,Y,f2(X,Y))  % surf plots surfaces in 3D 
% 
%% ridiculous amount of variables
% 
f3 = @(x,y,z,a,b,c) cos(x).*sin(y) + x.*log(z)...
                     - a./b + c.*cos(z);
f3(pts,pts,pts,pts,pts,pts)  % if I'm going to use 
                             % the function this way, I only need one
                             % variable
%% how many inputs?
% 
f3(1,2,3,4,5,6)   % this works
f3([1,2,3,4,5,6]) % this does not....why?
% 
%% what is happening here?
% 
f4 = @(x) x(1).^2 - x(2).^2.*x(3);  
% 
%% does it take scalars? vectors? rows or columns? matrices?
% 
%f4(5) 
%f4([1 2])
f4([1 2 3])
f4([1;2;3])
%f4(1,2,3)
f4([1 2 3 4 5 6 7 8 9])
%f4(1:9)
%f4([4 5 6; 7 8 9])    % what values does it use? 
% 
% %% same function with three input variables instead of one
% 
f5 = @(x,y,z) x.^2 - y.^2.*z;
f5(1,2,3)
% 
% %%  vector functions 
% 
f6 = @(x) [3*x.^2; 2*x-2; sin(x)];
f6(pts)

plot(pts,f6(pts))
xlabel('x'); ylabel('y')
title('three curves')
legend('3x^2','2x-2','sin(x)','Location','SouthEast')
