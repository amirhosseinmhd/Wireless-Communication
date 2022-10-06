function [b] = bit_gen(N, k)
rng(1)%%%%%%%%%%%%
b = randsrc(N, k, [0 1; 0.5 0.5]);
end

