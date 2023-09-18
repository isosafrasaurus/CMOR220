%Pierce Zhang, CMOR220, Fall 2023, Project 3: Finding real Roots
%week4_cooldrive.m
%This function implements Bisection and Newton's method to find the first root of the bar cooling problem for varying lengths and tolerances
%List of input: none
%List of output: none
%Last modified: September 11, 2023

function week4_cooldrive

    %% Part 1 (Bisection method for variable lenght)
    %intialize parameters
    a = 0.0001; %fixed left endpoint
    tol = 0.01; %fixed tolerance for figure 1
    L = 1:0.1:4; %variable length for figure 1
    x = zeros(1,length(L)); %preallocate x
    
    for cnt = 1:length(L)
        %correct the b (right endpoint) value so that there is only one root between the start and endpoint
        [x(cnt), ~] = debis(a,b,tol,L(cnt)); %run the debis function and compute row vector x
    end
    
    plot(L,x.^2,'x-') %plot cooling rate vs length
    %label your plot
    
    %% Part 2 (Compare Bisection and Newton's methods for variable tolerance)
    %open a new figure window for figure 2
    
    %initialize parameters
    L = 1; %fixed length for figure 2
    a = 0.1; %fixed left endpoint for figure 2
    b = 3; %fixed right endpoint for figure 2
    x = (a+b)/2; %midpoint - initial guess
    tol = 10.^-(1:8); %for the tolerance scale, create a vector of magnitudes
    iterB = zeros(1,length(tol)); %preallocate the iteration vector for bisection method
    iterN = zeros(1,length(tol)); %preallocate the iteration vector for Newton's method
    for i=1:length(tol)
       % [~,...] = debis(...) apply bisection method for specific tol
       % [~,...] = denewt(...) apply Newton's method for specific tol
    end
    
    semilogx(tol,iterB,'bx-')
    % include Newton's iterations vs Bisection
    
   % compute number of iterations required by Bisection (iterB) and Newton's methods (iterN) 
   % plot both iterB and iterN in one plot  
   % label your plot
end

function [x,iter] = denewt(x,tol,L)
%This function implements Newton's method and returns the root x and the
%number of iterations required to reach desirable tolerance
%List of input: x - an initial guess 
%               tol - tolerance 
%               L - length of the bar
%List of output: x - the root of the cooling function
%                iter - the number of iterations required to reach desirable tolerance using Newton's method


  %Add a check to make sure the derivative at the starting point is nonzero

  %Apply Newton's method


    
end

function [x, iter] = debis(a,b,tol,L)
%This function implements Bisection method and returns the root x and the
%number of iterations required to reach desirable tolerance
%List of input: a - fixed left endpoint
%               b - fixed right endpoint
%               tol - tolerance
%               L - length of the bar
%List of output: x - the root of the cooling function
%                iter - the number of iterations required to reach that desirable tolerance using bisection method

  %Add a check to make sure the root exists
  %Apply Bisection method

end

function val = coolfun(x,L)
%The function is f(x) = sin(xL) + xcos(xL)
%List of input: x - a parameter of the bar cooling model
%               L - length of the bar
%List of output: val - the current value of the bar cooling function

%Define the bar cooling function
end

function val = coolfundx(x,L)
%This function models the derivative w.r.t. x of the bar cooling problem
%List of input: x - a parameter of the bar cooling model
%               L - length of the bar
%List of output: val - the current value of the derivative of the bar cooling function

%Derivative of the bar cooling function
end