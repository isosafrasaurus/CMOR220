% Pierce Zhang, CMOR220, FALL 2023, Compentency nested loop and data
% comparison
% competency_nested_loops_data_comp.m
% Answers to nested loop and data comparison competency
% Last modified: October 21, 2023

% Driver to answer the questions
function competency_nested_loops_data_comp
    % Problem 1
    M = cone_matrix()

    clear;

    % Problem 2
    xs = -5:0.01:5;
    data = zeros(2, length(xs), 4);
    data(:,:,1) = Preprocess([xs;7.*sin(xs)]); %y1
    data(:,:,2) = Preprocess([xs;0.5.*exp(xs - 2)]); %y2
    data(:,:,3) = Preprocess([xs;4.*sqrt(xs + 5)]); %y3
    data(:,:,4) = Preprocess([xs;0.5 .* log(xs+6)]); %y4

    PWCM = zeros(4,4);
    for n1=1:4
        for n2=1:4
            PWCM(n1,n2)=DistanceCalc(data(:,:,n1),data(:,:,n2));
        end
    end
    figure();
    imagesc(PWCM); hold on;
    xlabel("First curve, y_{n2}"); ylabel("Second curve, y_{n2}");
    title("Problem 2: PWCM for curves");
    colorbar;
    
    clear;

    % Problem 3
    k = 5:5:25;
    tol = 10.^-(1:8);

    PWCM = zeros(length(k),length(tol));
    for k_i = 1:length(k)
        for tol_i = 1:length(tol)
            [~,PWCM(k_i,tol_i)] = Bisection(0,50,k(k_i),tol(tol_i));
        end
    end
    figure();
    imagesc(PWCM); hold on;
    xlabel("k, index"); ylabel("tol, index");
    title("Problem 3: PWCM for Bisection Function varying k,tol");
    colorbar;

    clear;

    % Problem 4
    b = 20:2:50;
    tol = 10.^-(1:12);

    PWCM = zeros(length(b),length(tol));
    for b_i = 1:length(b)
        for tol_i = 1:length(tol)
            [~,PWCM(b_i,tol_i)] = Bisection(0,b(b_i),10,tol(tol_i));
        end
    end
    figure();
    imagesc(PWCM); hold on;
    xlabel("b, index"); ylabel("tol, index");
    title("Problem 4: PWCM for Bisection Function varying b,tol");
    colorbar;
end

% Inputs: a, start interval, b, end interval, k, inverse multiplier for the
% leading term in the sample function, tol, tolerance
% Outputs: x, root value estimations, iter, number of iterations
function [x,iter]=Bisection(a,b,k,tol)
    % Function to execute the bisection method according to the
    % specification of problem 3.
    Eq=@(x,k) -(1/k)*x^2+x+10;
    x=(a+b)/2;
    iter=0;
    while abs(Eq(x,k))>tol
        if Eq(a,k)*Eq(x,k)>0
            a=x;
        else
            b=x;
        end
        x=(a+b)/2;
        iter=iter+1;
    end
end

% Inputs: C, vector of points to center and norm in preprocess
% Outputs: pC, preprocessed C
function [pC]=Preprocess(C)
    % Implements preprocessing for accurate and meaningful distance
    % measurement of curves. From RLA session.
    [M]=mean(C,2);
    centerC=C-M;
    pC=centerC/norm(centerC,'fro');
end

% Inputs: Cv1, vector of first curve, Cv2, vector of second curve
% Outputs: D, distance between the curves as specified in doc
function [D]=DistanceCalc(Cv1,Cv2)
    % Implements rotation elimination and distance calculation for
    % measurement of curves. From RLA session.
    A=Cv1*Cv2';
    [U,~,V]=svd(A);
    if det(U*V')>0
        D=norm(Cv1-(U*V')*Cv2,'fro');
    else
        D=norm(Cv1-(U*[1 0; 0 -1]*V')*Cv2,'fro');
    end
end

% Inputs: none
% Outputs: M, 3D matrix containing cone volume for specified ranges of
% a,b,h values.
function [M] = cone_matrix()
    M = zeros(4,4,4); a = 1:4; b = 5:8; h = 9:12;
    for a_i = 1:4
        for b_i = 1:4
            for h_i = 1:4
                M(a_i,b_i,h_i) = a(a_i)*b(b_i)*h(h_i)*pi/3;
            end
        end
    end
end