function [det_sym] = min_dist_detector(rx_sym, constellation)
[~, det_sym] = min(abs(rx_sym - constellation.').');
det_sym = det_sym.';