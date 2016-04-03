n = 50;
[A, b, x0] = get_A_b( n );
xexact = A\b;
w = 1.1;
TOL = 1 * 10^-15;

N = 100;

[x, residue] = SOR(n, A, b, x0, w, TOL, N);
disp(norm(x - xexact, inf));


[A, b, x0] = get_A_b( n );
[x] = CG( A, b, x0, N, TOL); 
disp(norm(x - xexact, inf));