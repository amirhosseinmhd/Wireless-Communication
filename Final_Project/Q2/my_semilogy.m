function [] = my_semilogy(x, y, title_name, x_label, y_label, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure;
semilogy(x , y)
grid on;
xlabel(x_label);
ylabel(y_label);
title(title_name);

if nargin>5
    varargin{:};
end

