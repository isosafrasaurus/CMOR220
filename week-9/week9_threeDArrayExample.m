% Anastasiya Protasov, CMOR220, Fall 2023, 3D matrices
% week9_threeDArrayExample.m
% Script to show examples of constructing and using 3-D arrays
% list of input: none
% list of output: none
% Last Modified: October 18, 2023

A = zeros(2,6,10);   % preallocate the array
A(:,:,1) = ones(2,6);  % the first page is an all-ones matrix
for i = 2:10
    A(:,:,i) = i*A(:,:,1);  % what are we doing?
end

A
 
%% what can we do with A? 

A(:)     % turn A into a single column vector of length i*j*k
pause
%A(:,:)   % turn A into a single matrix of dimensions i X j*k 

%% Arithmetic

B = ones(size(A));
clc
A+B                   % can do
pause
clc
A-B                 % can do
pause
clc
A.*B                % can do
pause 
clc
A./(2*B)              % can do
% pause
% A*B                   % not defined for 3d arrays
% pause
% A/B                   % not defined for 3d arrays
% pause 
% A\B                   % not defined for 3d arrays
pause
A.^2                  % can do
 
%% what can we not do? 

% A(:,10,10) = []   % ERROR
% A'                % ERROR
% norm(A)           % ERROR
% plot(A)           % ERROR 
% inv(A)            % ERROR

%% can we still transpose each page?

permute(A, [2,1,3])     % transposes each page
pause
permute(A,[3,1,2])      % restructure into a different array
pause
permute(A,[2,3,1])      % restructure into a different array

%% we can do some of the others using loops

for i = 1:10
    norm(A(:,:,i))   % take the matrix 2 norm of each page
end

pause
clc
for i = 1:10 
    norm(A(:,:,i), 'fro')  % take the frobenius norm of each page
end

pause
clc
for i = 1:10
    subplot(2,5,i)
    plot(A(:,:,i))    % plot each page
end

pause 
clc
for i = 1:10
    for j = 1:10
        A(:,:,i)*A(:,:,j)'   % computing 100 2X2 matrices
    end
end


