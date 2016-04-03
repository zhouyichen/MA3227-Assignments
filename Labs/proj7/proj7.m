p = @(x) x;
q = @(x) 4;
f = @(x) (4 * x.^2 - 8 * x + 1);
n = 20;
results  = RayleighRitz( p, q, f, n );