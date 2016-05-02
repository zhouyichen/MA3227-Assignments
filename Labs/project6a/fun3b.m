function f = fun3b(t,y)
[n,m]=size(y); f=zeros(n,m);
f(1) = y(2);
f(2) = -2/t*y(2)+2/(t^2)*y(1);
return
