function [xp_new, log] = ao_iter_once(P, xp, PT, opts)
%AO_ITER_ONCE  One outer AO iteration: P6 -> P8 -> P6 (fixed PT).
%   opts fields (all optional):
%     .solver   : 'sedumi' or 'sdpt3' (default 'sedumi')
%     .rho_pre  : Frobenius penalty for first P6 (default 10)
%     .rho_post : Frobenius penalty for second P6 (default 10 or 20)
%     .rho_i    : DC slack penalty (default 1)
%     .dc_iter  : DC iterations per inner solve (default 2)
%
% Outputs:
%   xp_new : updated positions after one sweep (P8)
%   log    : struct with progress signals

    if nargin<4, opts = struct(); end
    if ~isfield(opts,'solver'),   opts.solver   = 'sedumi'; end
    if ~isfield(opts,'rho_pre'),  opts.rho_pre  = 10; end
    if ~isfield(opts,'rho_post'), opts.rho_post = 10; end
    if ~isfield(opts,'rho_i'),    opts.rho_i    = 1; end
    if ~isfield(opts,'dc_iter'),  opts.dc_iter  = 2; end

    % --- Build geometry & primitives at current positions
    S  = pass.geom_channel(P, xp);
    B  = pass.beamformers(P, S);
    D  = pass.sdp_primitives(P, B, PT);

    % --- First inner solve (P6) gives targets
    ip1 = struct('rho',opts.rho_pre,'rho_i',opts.rho_i,'max_dc_iter',opts.dc_iter,...
                 'solver',opts.solver,'verbose',false);
    [~, ~, wmax, vmax, hist1, feas1] = pass.inner_sdp_dc_solve(P, D, ip1);

    % Align global phase (not strictly necessary, but helps per-element match)
    phi = -angle(wmax(1));
    wtilde = wmax * exp(1j*phi);
    vtilde = vmax * exp(1j*phi);

    % --- P8: one sweep of element-wise updates with fixed targets
    [xp_new, st] = pass.update_positions_sweep(P, xp, wtilde, vtilde);

    % --- Rebuild at new positions & solve inner again (P6)
    S2 = pass.geom_channel(P, xp_new);
    B2 = pass.beamformers(P, S2);
    D2 = pass.sdp_primitives(P, B2, PT);

    ip2 = struct('rho',opts.rho_post,'rho_i',opts.rho_i,'max_dc_iter',opts.dc_iter,...
                 'solver',opts.solver,'verbose',false);
    [Wtilde2, Vtilde2, wmax2, vmax2, hist2, feas2] = pass.inner_sdp_dc_solve(P, D2, ip2);

    % --- Progress signals (use raw sensing terms without penalties)
    raw_sense_before = real(trace(D.E  * D.Vx));     % at xp (model)
    raw_sense_after  = real(trace(D2.E * D2.Vx));    % at xp_new (model)
    raw_sense_opt    = real(trace(D2.E * Vtilde2));  % optimized inner at xp_new

    log = struct();
    log.dx_L2     = norm(xp_new - xp);
    log.min_space = min(diff(xp_new));
    log.box_ok    = (xp_new(1) >= -P.L/2 - 1e-9) && (xp_new(end) <= P.L/2 + 1e-9);

    log.raw_sense_before = raw_sense_before;
    log.raw_sense_after  = raw_sense_after;
    log.raw_sense_opt    = raw_sense_opt;

    log.feas1 = feas1;   log.hist1 = hist1;  % first inner
    log.feas2 = feas2;   log.hist2 = hist2;  % second inner
    log.Wtilde2 = Wtilde2; log.Vtilde2 = Vtilde2; log.wmax2 = wmax2; log.vmax2 = vmax2;
end