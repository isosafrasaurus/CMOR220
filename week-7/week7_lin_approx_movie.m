% Anastasiya Protasov, CMOR220, Fall 2023, movie example
% week7_lin_approx_movie.m
% Script to show and example of creating an animation in MATLAB using
% linear approximation
% list of input: none
% list of output: none
% Last Modified: October 4, 2023

function plot_lin_approx_movie

R = @(x) sin(5*x);
deriv = @(x) 5*cos(5*x);

delta = 1/20;
dom = 0:delta:3;

approx = zeros(1,length(dom));
approx(1) = R(0);

for i = 2:length(dom)
    approx(i) = approx(i-1) + delta*deriv(dom(i-1));
end

myvideo = VideoWriter('LinearApproxMovie.avi');
open(myvideo);
for ii = 1:length(dom)
    plot(0:0.001:3,R(0:0.001:3),'color','b'); hold on
    plot(dom(1:ii),approx(1:ii),'-o', 'color','r')
 %  plot(dom(ii),R(dom(ii)),'o');
    xlabel('time')
    ylabel('Population')         
    legend('R(t)', 'linear approximation')
    title(['time: ' num2str(dom(ii))]);
    currFrame  = getframe(gcf);
    writeVideo(myvideo,currFrame);
end
close(myvideo);
end