function [x,y]=FD_nonlinearBVP(xspan, alpha, beta, Nh, TOL, Mmax)
% 
% [x,y]=FD_nonlinearBVP(xspan,alpha,beta, Nh,Mmax) solve
% the differential equation
%             y'' = f(x, y, y')
% on the interval xspan(1) < x < xspan(2) with boundary conditions
%         y(xspan(1)) =alpha, y(xspan(2))=beta
% using the finite difference method on an equispaced grid of Nh+1 intervals.
%
% The functions f(x,y,u), and its partitial derivatives 
% fy(x,y,u), fu(x,y,u) are defined at the end of this file
%

a=xspan(1); b=xspan(2);
h=(b-a)/(Nh+1);  % divide the interval into Nh+1 sub-intervals

x=linspace(a, b, Nh+2); 
%x(1) and x(Nh+2) are boundary points; x(2:Nh+1) are interior points

% initial approximation
y = alpha + (beta-alpha)/(b-a)*(x'-a);

k=1;
while k <= Mmax
  u(2:Nh+1)= (y(3:Nh+2) - y(1:Nh))/(2*h); % first derivative of y at interior grid points
  
  J =zeros(Nh, Nh);
  d = zeros(Nh,1);
% first row of the Jacobi matrix and RHS vector 
  J(1,1) = 2 + h^2*fy(x(2),y(2),u(2));  
  J(1,2) = -1 + h/2 *fu(x(2), y(2), u(2));
  d(1) = -(2*y(2) - y(3) -alpha + h^2*f(x(2),y(2),u(2)));
  
% row 2 to row Nh-1
  for i=2:Nh-1
    J(i,i) = 2 + h^2*fy(x(i+1), y(i+1), u(i+1));
	tmp=fu(x(i+1), y(i+1), u(i+1));
	J(i,i+1)=-1+ h/2*tmp;
	J(i,i-1)=-1- h/2*tmp;
	d(i)=-(2*y(i+1)-y(i+2)-y(i)+ h^2*f(x(i+1), y(i+1), u(i+1)));
  end

% last row 
  J(Nh,Nh) = 2 + h^2*fy(x(Nh+1),y(Nh+1),u(Nh+1));
  J(Nh, Nh-1)= -1 - h/2*fu(x(Nh+1), y(Nh+1), u(Nh+1));
  d(Nh) = -(2*y(Nh+1)-y(Nh)-beta+ h^2*f(x(Nh+1), y(Nh+1), u(Nh+1)));

% solve the linear system
  v=J\d;
  
% update the approximation
  y(2:Nh+1) = y(2:Nh+1) + v;
  
% check convergence
  if norm(v, 2) <= TOL
    %fprintf('converged in %4d steps\n', k);
	return
  end
  
  k=k+1;
end

fprintf('failed to converge; max number of iterations exceeded\n');
return

function f = f(x, y, u)
f= u * cos(x) - y * log(y);
return

function f = fy(x, y, u)
f= - 1 - log(y);
return

function f = fu(x, y, u)
f= cos(x);
return
