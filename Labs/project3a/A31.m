g1 = @(x)((cos(x(2) * x(3)) + 0.5)/3);
g2 = @(x)((sqrt(x(1) ^ 2 + 0.3125)/25 - 0.03));
g3 = @(x)(exp(-x(1)*x(2))/(-20) - (10 * pi -3)/60);

G = {g1, g2, g3};

x0 = [0;0;0];
TOL = 10^-5;
N = 20;
x = FPI( x0, TOL, N, G);