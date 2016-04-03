function [ x ] = FPI( x0, TOL, N, G)
%FPI Summary of this function goes here
%   Detailed explanation goes here
    k = 1;
    n=size(x0);
    x = x0;
    while k <= N
        for i = 1:n
            x(i) = G{i}(x0);
        end
        if norm(x - x0, inf) <= TOL
            fprintf('\n FP converged after %5i iterations\n', k);
            return
        end
        x0 = x;
        k = k + 1;
    end
    fprintf('\n FP Maximum number (N =%5i) of iterations exceeded\n', N);
    return;
end

