function [det_sym_idx, rx_sym] = pulse_demodulation(rx_smpl, modulation, M, fs, smpl_per_symbl, pulse_name , rx_mode, varargin)

if modulation == "fsk"
    if nargin > 7
        mod_det_opt = varargin{1};
        if mod_det_opt == "coherent"
            coeff = 0.5;
        elseif mod_det_opt == "noncoherent"
            coeff = 1;
        end
    end
    rx_sym = [];
    for m = 0:M-1
        p = 1/sqrt(smpl_per_symbl) * exp(1i * 2*pi * ...
            m .* (0:smpl_per_symbl-1).' * 1/smpl_per_symbl * coeff);
        [rx_sym] = [rx_sym corr_match(rx_smpl, p, smpl_per_symbl, rx_mode)];
    end
    if nargin > 7
        mod_det_opt = varargin{1};
        if mod_det_opt == "coherent"
            [~, det_sym_idx] = max(real(rx_sym), [], 2);
        elseif mod_det_opt == "noncoherent"
            [~, det_sym_idx] = max(abs(rx_sym), [], 2);
        end
    end
    
    
else    
    pulse = pulse_shape(pulse_name, fs, smpl_per_symbl);
    [rx_sym] = corr_match(rx_smpl, pulse, smpl_per_symbl, rx_mode);
    [cons, ~] = constellation(M, modulation);
    det_sym_idx = min_dist_detector(rx_sym, cons);
end
