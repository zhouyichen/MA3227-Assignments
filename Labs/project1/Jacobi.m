% The Jacobi method for linear systems

% Input: the coefficient matrix A; RHS vector b; initial guess x0
%        tolerance TOL; max number of iterations N

% Output: approximate solution x, or an error message 

function [x, residue, error] = Jacobi(A, b, x0, TOL, N, xexact)

k=1;  % set the counter for iterations

n=size(A,1); % size of the linear system

while k <= N

 for i=1:n
    x(i,1)= ( -A(i,1:i-1)*x0(1:i-1) - A(i,i+1:n)*x0(i+1:n) +b(i))/A(i,i);
 end
 
 error(k)=norm(x-xexact, Inf);
 residue(k) = norm(b-A*x, Inf)/norm(b, Inf);  % residue vector
 if  residue(k) < TOL
   fprintf('\n converged after %5i iterations\n', k);
   return;
 end
 
 x0=x;
 k=k+1;
end

fprintf('\n Maximum number (N =%5i) of iterations exceeded\n', N)
