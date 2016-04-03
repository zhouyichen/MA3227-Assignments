function [ x ] = newton( x0, TOL, N, F, J )
%NEWTON Summary of this function goes here
%   Detailed explanation goes here
    k = 1;
    n=size(x0);
    x = x0;
    Fx = zeros(n, 1);
    Jx = zeros(n);
    while k <= N
        for i = 1:n
            Fx(i) = F{i}(x0);
            for j = 1:n
                Jx(i, j) = J{i}{j}(x0);
            end
        end
        y = Jx\(-Fx);
        x = x + y;
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

