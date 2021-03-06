f1 = @(x)(4 * x(1) - x(2) + x(3) - x(1)*x(4));
f2 = @(x)(-x(1) + 3 * x(2) - 2 * x(3) - x(2) * x(4));
f3 = @(x)(x(1) - 2 * x(2) + 3 * x(3) - x(3) * x(4));
f4 = @(x)(x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 - 1);
F = {f1, f2, f3, f4};

j11 = @(x)(4 - x(4));
j12 = @(x)(-1);
j13 = @(x)(1);
j14 = @(x)(2 * x(1));

j21 = @(x)(-1);
j22 = @(x)(3 - x(4));
j23 = @(x)(-2);
j24 = @(x)(2 * x(2));

j31 = @(x)(1);
j32 = @(x)(-2);
j33 = @(x)(3 - x(4));
j34 = @(x)(2 * x(3));

j41 = @(x)(-x(1));
j42 = @(x)(-x(2));
j43 = @(x)(-x(3));
j44 = @(x)(0);

j1 = {j11, j21, j31, j41};
j2 = {j12, j22, j32, j42};
j3 = {j13, j23, j33, j43};
j4 = {j14, j24, j34, j44};

J = {j1, j2, j3, j4};

x0 = [0;1;1;1];
TOL = 10^-5;
N = 20;
x = newton( x0, TOL, N, F, J);