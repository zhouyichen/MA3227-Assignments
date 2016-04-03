xspan = [1, 2];
f_exact = @(x)(4/x - 2/(x^2) + log(x) - 1.5);

N = 100;
step_sizes = zeros(N,1);
errors = zeros(N,1);
for s = 1:N
    n = 20 * s;
    step_size = (xspan(2) - xspan(1))/n;
    [x,y]=FD_linearBVP(xspan, 1/2, log(2), n);
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