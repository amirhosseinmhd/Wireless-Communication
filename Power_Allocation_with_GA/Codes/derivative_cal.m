function [dfR_dp] = derivative_cal(p, N0, G)

% k diffrent persons, i diffrent paths
diag_G = diag(G).';

% we have to sum over i other person(rather k) so we put i persons on
% columns.
cons = N0 + p*G' - p.*diag_G;
gamma_k = p.*diag_G./cons;
d_snr_equall_k = (diag_G./cons)./((1+ gamma_k).*(log10(1 + gamma_k)));


% p_i Gii * G_ik for i != k ; we repeat p_i G_i in order to have k 
% element vecctor
% rep_gp = repmat(, 50, 1);
d_snr_notequall_k_nominator = ((diag_G.*p)' .* G); % it will be summed over columns
% so it will be 1 50


% 
% rep_gp = repmat(diag_G.*p, 50, 1);
% d_snr_notequall_k_nominator = + (rep_gp * G).'; % it will be summed over columns
% % so it will be 1 50




% sum(p .* G, 2) , size = 50 1 
d_snr_notequall_k_denominator = (N0' + sum(p .* G, 2) - (p .* diag_G)').^2;
gamma_i = gamma_k';
R_i = log10(1 + gamma_i);

% d_snr_notequall_k = sum(d_snr_notequall_k_nominator, 1)./...
%     sum((d_snr_notequall_k_denominator)./(1 + gamma_i)./R_i);
x = d_snr_notequall_k_nominator./ ...
    d_snr_notequall_k_denominator./ (1 + gamma_i)./ R_i;
d_snr_notequall_k = sum(x,1 );



residual = (p.*diag_G.^2./(cons.^2))./(1+ gamma_k)./(log10(1 + gamma_k));

dfR_dp = d_snr_equall_k-(d_snr_notequall_k - residual); 
end

