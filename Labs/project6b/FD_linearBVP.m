function [x,y]=FD_linearBVP(xspan, alpha, beta,Nh)
% 
% [x,y]=FD_linearBVP(xspan,alpha,beta, Nh) solve
% the differential equation
%             y'' =p(x)y'+ q(x) y +r(x)
% on the interval xspan(1) < x < xspan(2) with boundary conditions
%         y(xspan(1)) =alpha, y(xspan(2))=beta
% using the fintie difference method on an equispaced grid of Nh+1 intervals.
%
% The functions p(x), q(x) and r(x) are defined at the end of this file
%

h=(xspan(2)-xspan(1))/(Nh+1);  % divide the interval into Nh+1 sub-intervals

x=linspace(xspan(1), xspan(2), Nh+2); 
%x(1) and x(Nh+2) are boundary points; x(2:Nh+1) are interior points

% evaluate p, q, r at the interior grid points
% 
p = fp(x(2:Nh+1));
q = fq(x(2:Nh+1));
r = fr(x(2:Nh+1));

A = zeros(Nh, Nh);  % coefficient matrix
for i=1:Nh
  A(i,i) = 2 + h^2*q(i);
end
for i=1:Nh-1
  A(i,i+1) = -1 + 1/2*h*p(i);
  A(i+1,i) = -1 - 1/2*h*p(i+1);
end

b=zeros(Nh,1); % creat a column vector 
               % for the right-hand side of the linear system
b(1) = -h^2 * r(1) + (1+1/2*h*p(1))*alpha;
b(Nh) = -h^2 * r(Nh) + (1-1/2*h*p(Nh))*beta;
for i=2:Nh-1
  b(i) = -h^2*r(i);
end

y=zeros(Nh+2,1);
y(1) = alpha;
y(2:Nh+1) = A\b;
y(Nh+2)=beta;

return

function f = fp(x)
f=-4./x;
return

function f = fq(x)
f=-2./(x.^2);
return

function f = fr(x)
f=2*(log(x))./(x.^2);
return
