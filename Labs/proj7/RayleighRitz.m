function [ results ] = RayleighRitz( p, q, f, n )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    h = 1/n;
    Q = zeros(6, n + 1);
    for i = 1:n+1
        xim1 = h * (i-1);
        xi = h * i;
        xip1 = h * (i+1);
        if (i < n)
            Q(2, i) = n^2 * integral(@(x)(((x - xim1).^2).*q(x)), xim1, xi);
            Q(3, i) = n^2 * integral(@(x) (((xip1 - x).^2).*q(x)), xi, xip1);
            Q(5, i) = n * integral(@(x) ((x - xim1).*f(x)), xim1, xi);
            Q(6, i) = n * integral(@(x) ((xip1 - x).*f(x)), xi, xip1);
            if (i < n - 1)
                Q(1, i) = n^2 * integral(@(x) ((xip1 - x).*(x - xim1).*q(x)), xi, xip1);
            end
        end
        Q(4, i) = n^2 * integral(@(x)(p(x)), xim1, xi);

    end
    A = zeros(n, n);
    B = zeros(n, 1);
    for i = 1:n
        A(i, i) = Q(4, i) + Q(4, i+1) + Q(2, i) + Q(3, i);
        if (i < n)
            A(i, i + 1) = - Q(4, i + 1) + Q(1, i);
        end
        if (i > 1)
            A(i, i - 1) = -Q(4, i) + Q(1, i - 1);
        end
        B(i) = Q(5, i) + Q(6, i);
    end
    C = linsolve(A, B);
    results = C;
    return
end
        

