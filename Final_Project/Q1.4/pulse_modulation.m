function [tx_smpl, cons] = pulse_modulation(sym_idx, modulation, M, fs, smpl_per_symbl, pulse_name , pulse_shape_mode, varargin)

if modulation == "fsk"
    if nargin > 7
        mod_det_opt = varargin{1};
        if mod_det_opt == "coherent"
            slm = 1/sqrt(smpl_per_symbl) * exp(1i * 2*pi * ...
                (0:M-1) .* (0:smpl_per_symbl-1).' * 1/smpl_per_symbl * 0.5);
        elseif mod_det_opt == "noncoherent"
            slm = 1/sqrt(smpl_per_symbl) * exp(1i * 2*pi * ...
                (0:M-1) .* (0:smpl_per_symbl-1).' * 1/smpl_per_symbl);
        end
    end
    tx_smpl = slm(:, sym_idx+1);
    tx_smpl = tx_smpl(:);
    cons = nan;
else    
    [cons, ~] = constellation(M, modulation);
    mod_sym = cons(sym_idx+1);
    
    pulse = pulse_shape(pulse_name, fs, smpl_per_symbl);
    if pulse_shape_mode == "kron"
        tx_smpl = kron(mod_sym, pulse);
    elseif pulse_shape_mode == "conv"
        s_zero_pad = upsample(mod_sym, smpl_per_symbl);
        s_zero_pad = s_zero_pad(1:end-smpl_per_symbl+1);
        tx_smpl = conv(s_zero_pad, pulse);
    end
end