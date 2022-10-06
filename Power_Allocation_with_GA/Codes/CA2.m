 clc
clear
rng(1)
close 
% for j == 4 change learning rate to 70 for convergence
for j = [1 2 3 5]
    filename = sprintf('a%d.mat', j);
    load(filename);
    [~, N] = size(P_max);
    a = 0;
    b = max(P_max);
    p = a + (b-a).*rand(1,N);
%     p = ones(1, 50);
    f_r = @(P) sum(log10(log10(1 +  P.*(diag(G).')./(N0 + P*G.' - P.*(diag(G).')))));
    f_prev = f_r(p);
    f = 0;
    list_fr = [];
    alpha = 10;
    i = 0;
    while abs(f - f_prev) > 10^-4
        update_p = derivative_cal(p, N0, G);
        p = p + alpha * update_p;
        p(p>P_max) = P_max(p>P_max);
        p(p<0) = 0.000000001;
        f_prev = f;
        f = f_r(p);
        list_fr = [list_fr, f];
        i = i + 1;

    end
    
    %%
%     %%%%%%%%%%%%%%% Part B.1 uncomment here
%     figure
%     plot(list_fr)
%     filename = sprintf('Optimum Point of a%d.mat', j);
%     title(filename);
%     xlabel("Number of Iterations")
%     ylabel("f(r)")
%     grid on
    %%
%     %%%%%%%%%%%%% Part B.2
%     stem(p)
%     figure
%     filename = sprintf('Stem for of a%d.mat P^o', j);
%     title(filename);
%     %%
%     %%%%%%%%%%%%% Part B.3 Uncomment here
%     figure
%     stem(P_max - p)
%     filename = sprintf('Stem for of a%d.mat P^m_i - P^o', j);
%     title(filename);
%%
    f_max = f_r(P_max);
    f_optimum = f_r(p);
    increase = f_optimum - f_max;
    saved_power = sum(P_max - p);
    fprintf("We have %.4f Imporvement for Optimum value for a%d.mat and also %.2f Power has been saved\n"...
    ,increase, j, saved_power);
    

end
