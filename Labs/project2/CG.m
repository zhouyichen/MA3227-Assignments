function [ x ] = CG( A, b, x0, N, TOL)
%CG Summary of this function goes here
%   Detailed explanation goes here
    r = b - A * x0;
    v = r;
    x = x0;
    p = (r' * r);
    treshold = TOL * norm(b);
    for i = 1:N
        if sqrt(p) <= treshold
            fprintf('\n CG converged after %5i iterations\n', i);
            return;
        end
        A_v = A * v;
        t = p / (v' * A_v);
        x = x + t * v;
        r = r - t * A_v;
        new_p = (r' * r);
        s = new_p / p;
        p = new_p;
        v = r + s * v;
    end
    fprintf('\n CG Maximum number (N =%5i) of iterations exceeded\n', N)
    return;
end

