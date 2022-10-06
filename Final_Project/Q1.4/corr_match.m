function [rx_sym] = corr_match(rx_smpl, p, smpl_per_symbl, rx_mode)
if rx_mode == "correlator"
    rx_corr = xcorr(rx_smpl, p);
    rx_sym = rx_corr(length(rx_smpl):smpl_per_symbl:end);
elseif rx_mode == "matched_filter"
    rx_conv = conv(rx_smpl, flip(conj(p)));
    rx_sym = rx_conv(smpl_per_symbl:smpl_per_symbl:end);
end
