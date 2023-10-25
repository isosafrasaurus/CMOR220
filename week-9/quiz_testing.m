C1 = [-1 -0.5 0 0.5 1;-1 -0.5 0 0.5 1];
C2 = [-1 -0.5 0 0.5 1;1 0.5 0 -0.5 -1];

A = C1*C2';
[U,~,V] = svd(A);
deter = det(U*V');
if deter > 0
    O = U*V'
else
    O = U * [1 0;0 -1] * V'
end
F = C1 - U*V'*C2;
norm(F,'fro')