function [ A, b, x0] = get_A_b( n )
%GET_A_B Summary of this function goes here
%   Detailed explanation goes here
    A = eye(n);
    for i = 1:n-1
        A(i, i+1) = -1/2;
        A(i+1, i) = -1/2;
    end
    b = zeros(n, 1);
    b(1) = 1/2;
    x0 = zeros(n, 1);
end

