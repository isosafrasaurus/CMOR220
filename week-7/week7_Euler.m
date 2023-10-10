% Anastasiya Protasov, CAAM 210, Fall 2023, Euler's implementation
% week7_EulerExample.m
% Script to apply Euler's method to find an approximation of the solution
% for dx/dt = 2*x*t, x(0) = 8
% list of input: none
% list of output: none
% Last Modified: October 2, 2023


function Eulers
t=0:.02:4;
%x=zeros(length(t),1);
x(1)=8;
%delta = (t(end)-t(1))/length(t);
delta=0.01;
for k=1:length(t)-1 % dont put in first cell bro
    x(k+1)=x(k)+delta*(x(k)*(x(k)+t(k))); % do the euler thing
end
plot(t,x,'b-');
hold on
plot(t,8*exp(t.^2),'r-');
legend('approximation','solution')
end

