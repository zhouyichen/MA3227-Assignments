function [ A, b, x0] = get_A_b( n )
%GET_A_B Summary of this function goes here
%   Detailed explanation goes here
    A = zeros(n);
    for i = 1:n
        if i < n
            A(i, i+1) = -1;
            A(i+1, i) = -1;
        end
        A(i, i) = 2 * i;
    end
    b = 1.5 * (1:n) - 6;
    b = b';
    x0 = zeros(n, 1);
end

