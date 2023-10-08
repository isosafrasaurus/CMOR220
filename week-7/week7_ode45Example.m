% Anastasiya Protasov, CMOR220, Fall 2023, ode45 example
% week7_ode45Example.m
% Script to illustrate how to use ode45
% list of input: none
% list of output: none
% Last Modified: October 4, 2023


f = @(t,R) R*(1 - R/100);    % logistic growth model with r=1, K=100
[t,R] = ode45(f,[0,3],10);     % solve the ODE
hold on
plot(t,R,'-o','linewidth',2);  % plot the solution from ode45
plot(0:3, [10,19,34.39,56.953279],'-d','linewidth',2);
% actual solution
t1 = 0:0.01:3;
R1 = 1000./(10+90*exp(-t1));
plot(t1,R1,'-','linewidth',2);% solution using Euler's method with delta = 1
hold off
title('ode45 vs forward Euler')
xlabel('time')
ylabel('population')
legend('ode45 solution', 'forward Euler solution', 'actual solution', 'Location', 'NorthWest')
abs(R(end)-56.953279)

%% using options
options = odeset('RelTol', 1e-1);       % set the accuracy of your solution
[t2,R2] = ode45(f,[0,3],10,options);
figure 
plot(t,R,'-o', t2,R2,'-d','linewidth',2)
title('default vs custom tolerance')
xlabel('time')
ylabel('population')
legend('default tolerance', 'custom tolerance', 'Location', 'NorthWest')


%% longer timescale
% 
[t3, R3] = ode45(f,[0,10],10,options);
plot(t3,R3,'-o','linewidth',2)
title('Logistic Growth Solution for [0,10]')
xlabel('time')
ylabel('population')
% 

%% systems of ODE. 
options = odeset('RelTol', 1e-5); 
V = @(t,y) [y(2);(1-y(1).^2).*y(2) - y(1)];  % system representing a Van der Pol oscillator
[t1,y1] = ode45(V,[0,20],[2,0],options);     
plot(t1,y1,'linewidth',2)
title('Solution to the Van der Pol Equation')
xlabel('time')
ylabel('solution')

%% another system using an external function 
close all
clear all
lorenz  = @(t,y) [10*(y(2)-y(1));28*y(1)-y(2)-y(1).*y(3);y(1)*y(2)-8/3*y(3)];
[t2,y2] = ode45(lorenz,[0,70],[0,1,0]);
figure
plot3(y2(:,1),y2(:,2),y2(:,3))   % plot the curve defined by the three solutions
