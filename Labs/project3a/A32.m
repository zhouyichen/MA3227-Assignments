g1 = @(x)(sqrt(x(2)^2/5));
g2 = @(x)((sin(x(1)) + cos(x(2)))/4);

G = {g1, g2};

x0 = [0.25; 0.25];
TOL = 10^-5;
N = 20;
x = FPI( x0, TOL, N, G);