function [xp_new, log] = ao_iter_paper(P, xp, PT, rho, opts)
%AO_ITER_PAPER  One outer AO iteration (Algorithm 2) with a single rho.
%   - Run inner DC (P6) with rho to get targets
%   - One sweep of P8 (positions)
%   - Run inner DC (P6) again with the SAME rho
%   Returns xp_new and diagnostics incl. penalty-gap (chi2).
%
% opts fields (optional): solver ('sedumi'|'sdpt3'), dc_iter (default 2), rho_i (default 1)

    if nargin<5, opts = struct(); end
    if ~isfield(opts,'solver'),  opts.solver  = 'sedumi'; end
    if ~isfield(opts,'dc_iter'), opts.dc_iter = 2; end
    if ~isfield(opts,'rho_i'),   opts.rho_i   = 1; end

    % --- Build at current positions
    S  = pass.geom_channel(P, xp);
    B  = pass.beamformers(P, S);
    D  = pass.sdp_primitives(P, B, PT);

    % --- Inner DC (P6) with rho -> targets
    ip = struct('rho',rho,'rho_i',opts.rho_i,'max_dc_iter',opts.dc_iter,...
                'solver',opts.solver,'verbose',false);
    [~, ~, wmax, vmax, hist1, feas1] = pass.inner_sdp_dc_solve(P, D, ip);

    % Phase align targets (optional)
    phi = -angle(wmax(1));
    wtilde = wmax * exp(1j*phi);
    vtilde = vmax * exp(1j*phi);

    % --- P8 sweep with fixed targets
    [xp_new, st] = pass.update_positions_sweep(P, xp, wtilde, vtilde);

    % --- Rebuild & inner DC again with the SAME rho
    S2 = pass.geom_channel(P, xp_new);
    B2 = pass.beamformers(P, S2);
    D2 = pass.sdp_primitives(P, B2, PT);

    [Wtilde2, Vtilde2, wmax2, vmax2, hist2, feas2] = pass.inner_sdp_dc_solve(P, D2, ip);

    % --- Paper metrics: chi2 penalty-gap & sensing terms (model / optimized)
    chi2 = norm(Wtilde2 - D2.Wx, 'fro') + norm(Vtilde2 - D2.Vx, 'fro');
    raw_sense_before = real(trace(D.E  * D.Vx));    % model at xp
    raw_sense_after  = real(trace(D2.E * D2.Vx));   % model at xp_new
    raw_sense_opt    = real(trace(D2.E * Vtilde2)); % optimized inner at xp_new

    log = struct();
    log.dx_L2 = norm(xp_new - xp);
    log.min_space = min(diff(xp_new));
    log.box_ok = (xp_new(1) >= -P.L/2 - 1e-9) && (xp_new(end) <= P.L/2 + 1e-9);
    log.pen_gap = chi2;
    log.raw_sense_before = raw_sense_before;
    log.raw_sense_after  = raw_sense_after;
    log.raw_sense_opt    = raw_sense_opt;
    log.feas1 = feas1; log.hist1 = hist1;
    log.feas2 = feas2; log.hist2 = hist2;
    log.Wtilde2 = Wtilde2; log.Vtilde2 = Vtilde2; log.wmax2 = wmax2; log.vmax2 = vmax2;
end