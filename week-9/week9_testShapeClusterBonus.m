% Anastasiya Protasov, CMOR220, Fall 2023, use of linkage
% week9_testShapeClusterBonus.m
% Script to show an example of how to do clustering using linkage
% list of input: none
% list of output: none
% Last Modified: October 18, 2023

    X = randn(5, 2);
    figure(1);clf
    shapes = {'or', '*b', 'xk', 'vm', '^g'};
    for i = 1 : size(X, 1)
        plot(X(i, 1), X(i, 2), shapes{i}, 'markersize', 10);
        hold on
    end
    axis equal
    legend('1', '2', '3', '4', '5');
    pdist2(X, X)
    Z = linkage(pdist(X));
    c = cluster(Z, 'maxclust', 3)

