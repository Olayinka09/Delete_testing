function [xp_new, x_opt, f_opt] = update_position_1d(P, xp, n, wtilde_n, vtilde_n)
%UPDATE_POSITION_1D  Element-wise 1D position update for P8.
%   Minimizes |wtilde_n - w_n(x)|^2 + |vtilde_n - v_n(x)|^2
%   subject to box & spacing constraints on x in [-L/2, L/2].
%
% Inputs:
%   P         : config_pass()
%   xp        : (N x 1) current positions
%   n         : index (1..N) of the element to update
%   wtilde_n  : target complex scalar (nth entry of w~)
%   vtilde_n  : target complex scalar (nth entry of v~)
%
% Outputs:
%   xp_new : updated positions (same as xp except xp_new(n)=x_opt)
%   x_opt  : optimal x for element n
%   f_opt  : minimum objective value at x_opt

    N  = numel(xp);
    L  = P.L;
    dx = P.DELTAx;

    % Feasible interval from box & spacing constraints
    lb = -L/2;  ub =  L/2;
    if n > 1,   lb = max(lb, xp(n-1) + dx); end
    if n < N,   ub = min(ub, xp(n+1) - dx); end

    % If infeasible due to neighbors, keep position unchanged
    if lb > ub
        xp_new = xp; x_opt = xp(n); f_opt = inf; return;
    end

    % Precompute static geometry for user/target & feed
    psi0 = [0, 0, P.d];   % feed on the guide
    psi_c = [P.rc*cos(P.phic), P.rc*sin(P.phic), 0];
    psi_s = [P.rs*cos(P.phis), P.rs*sin(P.phis), 0];

    % Equal power per element
    alpha_n = (P.alpha_sum / P.N);

    % Objective handle: depends only on x for the nth component
    function f = obj_at(x)
        % nth element coordinate (x,0,d)
        p = [x, 0, P.d];

        % distances to user/target
        rc = sqrt(sum((p - psi_c).^2));
        rs = sqrt(sum((p - psi_s).^2));

        % in-waveguide phase from feed
        Lfeed = sqrt(sum((p - psi0).^2));
        theta = 2*pi * P.eta_eff * Lfeed / P.lambda;

        % n-th beamformer components
        w_n = exp(-1j*((2*pi/P.lambda)*rc + theta)) / (sqrt(alpha_n)*rc);
        v_n = exp(-1j*((2*pi/P.lambda)*rs + theta)) / (sqrt(alpha_n)*rs);

        % element-wise mismatch (P8)
        f = abs(wtilde_n - w_n)^2 + abs(vtilde_n - v_n)^2;
    end

    % Coarse-to-fine 1-D search (robust & simple)
    % Coarse grid
    Mcoarse = 301;
    xs = linspace(lb, ub, Mcoarse);
    fs = arrayfun(@obj_at, xs);
    [f_best, k] = min(fs);
    x_best = xs(k);

    % Local refine with small golden-section around x_best
    span = max(ub - lb, eps);
    a = max(lb, x_best - 0.05*span);
    b = min(ub, x_best + 0.05*span);
    phi = (sqrt(5)-1)/2;   % golden ratio conjugate
    c = b - phi*(b-a);
    d = a + phi*(b-a);
    fc = obj_at(c); fd = obj_at(d);
    for it = 1:40
        if fc < fd
            b = d;  d = c;  fd = fc;  c = b - phi*(b-a); fc = obj_at(c);
        else
            a = c;  c = d;  fc = fd;  d = a + phi*(b-a); fd = obj_at(d);
        end
        if (b - a) < 1e-4, break; end
    end
    if fc < fd, x_opt = c; f_opt = fc; else, x_opt = d; f_opt = fd; end

    % Apply update
    xp_new = xp;
    xp_new(n) = x_opt;
end