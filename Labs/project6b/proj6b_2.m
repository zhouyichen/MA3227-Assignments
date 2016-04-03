xspan = [0, pi/2];
f_exact = @(x)(exp(sin(x)));

N = 100;
step_sizes = zeros(N,1);
errors = zeros(N,1);
for s = 1:N
    n = 20 * s;
    step_size = (xspan(2) - xspan(1))/n;
    Mmax = 10;
    [x,y] = FD_nonlinearBVP(xspan, 1, exp(1), n, 10^-6, Mmax);
    step_sizes(s) = step_size;
    largest_error = 0;

    for i = 1:n
        x_c = x(i);
        error = abs(f_exact(x_c) - y(i));
        if (error > largest_error)
            largest_error = error;
        end
    end
    errors(s) = largest_error;
end
display(errors);
plot(step_sizes, errors);