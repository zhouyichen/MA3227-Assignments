function [t,u]=feuler(odefun, tspan, y0, Nh, varargin)
% 
% [T,Y]=feuler(odefun, tspan, y0, Nh) with tspan=[t0,tf] integrates the system
% of differential equations y'=f(t,y) from time t0 to tf with initial condition
% y0 using the forward Euler method on an equispaced grid of Nh intervals.
%
% Function odefun(t, y) must return a vector, whose elements hold the 
% evaluation of f(t,y), of the same dimension of y. 
%
% Each row in the solution array Y corresponds to a time returned in the 
% column vector T. 
%
% [T, Y]=feuler(odefun, tspan, y0, Nh, P1, P2,...) passes the additional
% parameters P1, P2, ... to the function odefun as odefun(t,y, P1, P2, ...).
%

h=(tspan(2)-tspan(1))/Nh;
y=y0(:); % always creates a column vector
w=y; 
u=y';

tt=linspace(tspan(1), tspan(2), Nh+1);

for t = tt(1:end-1)
  w=w+ h*feval(odefun, t, w, varargin{:});
  u = [u; w'];
end
t=tt;
return

