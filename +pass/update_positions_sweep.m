function [xp_new, stats] = update_positions_sweep(P, xp, wtilde, vtilde)
% One left-to-right sweep of element-wise updates (P8) with fixed targets.
    N = numel(xp);
    xp_new = xp;
    fvals  = zeros(N,1);
    for n = 1:N
        [xp_new, ~, fvals(n)] = pass.update_position_1d(P, xp_new, n, wtilde(n), vtilde(n));
    end
    stats = struct();
    stats.fvals = fvals;
    stats.min_spacing = min(diff(xp_new));
    stats.box_ok = (xp_new(1) >= -P.L/2 - 1e-9) && (xp_new(end) <= P.L/2 + 1e-9);
end