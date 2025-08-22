function [xp, out] = ao_run_paper(P, xp0, PT, rho0, c2, eps2, eps3, opts)
% AO_RUN_PAPER  Algorithm 2 (Penalty-based AO) exactly as in the paper.
%   repeat
%     repeat
%       update {w~, v~} by Algorithm 1 (P6 + DC)
%       update x^p via element-wise optimization (P8 sweep)
%     until objective value converges with accuracy eps2
%     rho <- c2 * rho
%   until ||W~-W(x^p)||_F + ||V~-V(x^p)||_F <= eps3   (measured in bar-space via tW+tV)
%
%   χ2 is measured in the paper's normalized scale:  χ2 := tW + tV  from Alg.1.

    if nargin < 8, opts = struct(); end
    if ~isfield(opts,'solver'),    opts.solver    = 'sedumi'; end
    if ~isfield(opts,'dc_iter'),   opts.dc_iter   = 3;  end         % DC iterations inside Alg.1
    if ~isfield(opts,'rho_i'),     opts.rho_i     = 1;  end         % initial DC slack penalty
    if ~isfield(opts,'c1bar'),     opts.c1bar     = 0.8; end        % Alg.1 Step 5 (shrink)
    if ~isfield(opts,'eps1'),      opts.eps1      = 1e-3; end       % Alg.1 Step 6 (slack stop)
    if ~isfield(opts,'max_outer'), opts.max_outer = 10; end
    if ~isfield(opts,'max_inner'), opts.max_inner = 20; end
    if ~isfield(opts,'verbose'),   opts.verbose   = true; end

    xp  = xp0(:);
    rho = rho0;

    % logs
    out.rho          = zeros(opts.max_outer,1);
    out.obj          = zeros(opts.max_outer,1);   % inner objective after inner-repeat
    out.chi2         = zeros(opts.max_outer,1);   % normalized penalty gap (tW+tV)
    out.dxL2         = zeros(opts.max_outer,1);   % last inner sweep Δx (L2)
    out.min_spacing  = zeros(opts.max_outer,1);
    out.box_ok       = false(opts.max_outer,1);

    for k = 1:opts.max_outer
        if opts.verbose
            fprintf('\n=== OUTER %d | rho=%.3g ===\n', k, rho);
        end

        % ----- INNER REPEAT (Alg.2 lines 3–6) -----
        obj_prev  = NaN;
        dx_last   = NaN;
        minsp_last= NaN;
        box_last  = false;
        chi2      = Inf;

        for j = 1:opts.max_inner
            % (A) Algorithm 1 at current xp -> w~, v~
            S  = pass.geom_channel(P, xp);
            B  = pass.beamformers(P, S);
            D  = pass.sdp_primitives(P, B, PT);

            ip = struct('rho',rho,'rho_i',opts.rho_i,'c1bar',opts.c1bar, ...
                        'eps1',opts.eps1,'max_dc_iter',opts.dc_iter, ...
                        'solver',opts.solver,'verbose',false);

            [~, ~, wmax, vmax, ~, ~] = pass.inner_sdp_dc_solve(P, D, ip);

            % phase-align targets
            phi    = -angle(wmax(1));
            wtilde = wmax * exp(1j*phi);
            vtilde = vmax * exp(1j*phi);

            % (B) P8: element-wise sweep to update positions
            [xp_new, st] = pass.update_positions_sweep(P, xp, wtilde, vtilde);

            % (C) Evaluate objective & normalized gap at new xp (Alg.1 once)
            S2 = pass.geom_channel(P, xp_new);
            B2 = pass.beamformers(P, S2);
            D2 = pass.sdp_primitives(P, B2, PT);
            [~, ~, ~, ~, histB, feasB] = pass.inner_sdp_dc_solve(P, D2, ip);

            obj  = histB.obj(end);                 % inner objective (scaled)
            chi2_abs  = histB.tW(end) + histB.tV(end);  % normalized penalty gap (paper scale)

            % relative change (skip test on first inner)
            if isfinite(obj_prev)
                relchg = abs(obj - obj_prev) / max(1, abs(obj_prev));
            else
                relchg = Inf;
            end

            if opts.verbose
                fprintf('  [inner %2d] obj=%.3e | Δrel=%.2e | χ2=%.3e (tW=%.3e,tV=%.3e) | QoS: %.2e>=%.2e\n', ...
                        j, obj, relchg, chi2, histB.tW(end), histB.tV(end), feasB.qos_lhs, feasB.qos_rhs);
            end

            % Δx, spacing, box BEFORE updating xp
            dx_last    = norm(xp_new - xp);
            minsp_last = st.min_spacing;
            box_last   = st.box_ok;

            % prepare next inner repeat
            xp       = xp_new;
            obj_prev = obj;

            % stop inner when objective converges to accuracy eps2
            if relchg < eps2
                break;
            end
        end

        % ----- record outer-iteration summary -----
        out.rho(k)         = rho;
        out.obj(k)         = obj_prev;
        out.chi2(k)        = chi2;
        out.dxL2(k)        = dx_last;
        out.min_spacing(k) = minsp_last;
        out.box_ok(k)      = box_last;

        if opts.verbose
            fprintf('  OUTER %d summary: obj=%.3e | χ2=%.3e | minΔ=%.3f | box=%d\n', ...
                    k, out.obj(k), out.chi2(k), out.min_spacing(k), out.box_ok(k));
        end

        % ----- OUTER STOP (Alg.2 line 8) -----
        if chi2 <= eps3
            if opts.verbose
                fprintf('Converged: penalty gap χ2 ≤ eps3 (%.3e ≤ %.3e)\n', chi2, eps3);
            end
            % truncate logs to k
            out.rho          = out.rho(1:k);
            out.obj          = out.obj(1:k);
            out.chi2         = out.chi2(1:k);
            out.dxL2         = out.dxL2(1:k);
            out.min_spacing  = out.min_spacing(1:k);
            out.box_ok       = out.box_ok(1:k);
            return;
        end

        % ----- STEP 8 (Alg.2 line 7): shrink rho -----
        rho = c2 * rho;
    end

    % truncate logs if max_outer reached
    lastk = find(out.rho ~= 0, 1, 'last');
    if ~isempty(lastk)
        out.rho         = out.rho(1:lastk);
        out.obj         = out.obj(1:lastk);
        out.chi2        = out.chi2(1:lastk);
        out.chi2_abs     = zeros(opts.max_outer,1);      % absolute gap (tW+tV)
        out.dxL2        = out.dxL2(1:lastk);
        out.min_spacing = out.min_spacing(1:lastk);
        out.box_ok      = out.box_ok(1:lastk);
    end
end