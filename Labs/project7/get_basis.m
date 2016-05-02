function [ basis ] = get_basis( n, h, i )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    function [y] = helper(x) 
        if (x <= h * (i -1))
            y = 0;
        elseif (x > h * (i -1) && x <= h * i)
            y = (x - h * (i -1))/h;
        elseif (x <= h * (i+1) && x > h * i)
            y = (x - h * (i -1))/h;
        else
            y = 0;
        end
        return;
    end;
    basis = helper;
    return;
end

