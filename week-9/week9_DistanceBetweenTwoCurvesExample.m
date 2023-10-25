% Anastasiya Protasov, CMOR220, Fall 2023, "similar" curves
% week9_DistanceBetweenTwoCurvesExample.m
% Script to show an example of how to compute a distance between two
% list of input: none
% list of output: none
% Last Modified: October 19, 2022
%% part 1: a rotation
D = [-2:0.05:2; (-2:0.05:2).^2]; % these points are on a parabola
plot(D(1,:), D(2,:), 'k.', 'MarkerSize',15)
pause
theta = pi/6; % angle of rotation
rot = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % rotation matrix
Drot = rot*D; % rotate the parabola
hold on
plot(Drot(1,:), Drot(2,:),'r.','MarkerSize', 15)
axis equal
hold off
%% part 2: removing the rotation
E = D*Drot';
[U,~,V] = svd(E);
det(U*V')
pause
F = D - U*V'*Drot;
norm(F,'fro')