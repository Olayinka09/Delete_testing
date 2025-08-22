function [Wtilde, Vtilde, wmax, vmax, hist, feas] = inner_sdp_dc_solve(P, D, opts)
%INNER_SDP_DC_SOLVE  Fixed-position inner solve (P6) with DC (rank-one) iterations.
%   Uses scaling of diagonal caps and normalization of E for numerical stability.
%   This implements the "inner layer" of Algorithm 1 for a fixed layout x^p.
%
% Inputs:
%   P    : config_pass() struct
%   D    : pass.sdp_primitives(P,B,PT) output (E, Wx, Vx, gammaQoS, diag caps, wmax0, vmax0, N)
%   opts : struct with optional fields:
%          - rho         (double, default 10)   : Frobenius mismatch penalty
%          - rho_i       (double, default 1)    : DC slack penalty
%          - max_dc_iter (int,    default 3)    : DC iterations
%          - tol_rel     (double, default 1e-3) : relative obj improvement tol
%          - solver      (char,   default 'sedumi') % 'sedumi' or 'sdpt3'
%          - verbose     (logical, default false)
%
% Outputs:
%   Wtilde,Vtilde : (N x N) Hermitian PSD in ORIGINAL scale (last DC iterate)
%   wmax,vmax     : principal eigenvectors of Wtilde,Vtilde (unit-norm)
%   hist          : struct with per-iter logs (obj, rank proxies, tW,tV,varpi,...)
%   feas          : struct with final feasibility diagnostics

    % ------------------ defaults & options ------------------
    if ~exist('opts','var'), opts = struct(); end
    def = struct('rho',10, 'rho_i',1, 'max_dc_iter',3, 'tol_rel',1e-3, ...
                 'solver','sedumi', 'verbose',false);
    fn = fieldnames(def);
    for k = 1:numel(fn)
        if ~isfield(opts,fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
    end

    % ------------------ scaling / normalization -------------
    N  = D.N;
    sW = D.diagW_cap(1);    % > 0
    sV = D.diagV_cap(1);
    if sW <= 0 || sV <= 0, error('Diagonal caps must be positive.'); end

    % QoS RHS in bar-space
    rhs_qos  = D.gammaQoS * P.sigma2 / sW;

    % Normalize E by its trace to reduce conditioning
    tE  = trace(D.E) + eps;
    En  = D.E / tE;
    rhs_qos_n = rhs_qos / tE;

    % Targets in bar-space (for Frobenius penalties)
    Wbar0 = D.Wx / sW;
    Vbar0 = D.Vx / sV;

    % Initial DC linearization refs (unit-norm)
    wref = D.wmax0 / max(norm(D.wmax0), eps);
    vref = D.vmax0 / max(norm(D.vmax0), eps);

    % ------------------ history containers ------------------
    hist.obj      = zeros(opts.max_dc_iter,1);
    hist.rankW    = zeros(opts.max_dc_iter,1);
    hist.rankV    = zeros(opts.max_dc_iter,1);
    hist.tW       = zeros(opts.max_dc_iter,1);
    hist.tV       = zeros(opts.max_dc_iter,1);
    hist.varpi1   = zeros(opts.max_dc_iter,1);
    hist.varpi2   = zeros(opts.max_dc_iter,1);
    hist.qos_lhs  = zeros(opts.max_dc_iter,1);
    hist.qos_rhs  = rhs_qos_n * ones(opts.max_dc_iter,1);

    % ------------------ DC iterations -----------------------
    obj_prev = -inf;
    Wbar = []; Vbar = [];
    for it = 1:opts.max_dc_iter
        % --- choose solver / precision ---
        cvx_clear
        switch lower(opts.solver)
            case 'sedumi'
                cvx_solver sedumi
            otherwise
                cvx_solver sdpt3
                cvx_solver_settings('rmdepconstr',1)
        end
        if opts.verbose, cvx_begin sdp; cvx_precision best
        else,            cvx_begin sdp quiet; cvx_precision best
        end

            % Variables (bar-space)
            variable Wbar(N,N) hermitian semidefinite
            variable Vbar(N,N) hermitian semidefinite
            variable tW nonnegative
            variable tV nonnegative
            variable varpi1 nonnegative
            variable varpi2 nonnegative

            % Objective: sensing - Frobenius penalties - DC slack penalties
            maximize( sV * trace(En * Vbar) ...
                      - (1/(2*opts.rho))  * (tW + tV) ...
                      - (1/(2*opts.rho_i)) * (varpi1 + varpi2) )

            subject to
                % Epigraphs (Frobenius norms)
                norm(Wbar - Wbar0, 'fro') <= tW
                norm(Vbar - Vbar0, 'fro') <= tV

                % Diagonal caps
                diag(Wbar) <= 1
                diag(Vbar) <= 1

                % QoS (normalized & scaled)
                trace(En * Wbar) >= rhs_qos_n

                % DC rank-one linearizations about (wref, vref)
                real( trace( Wbar * (eye(N) - (wref*wref')) ) ) <= varpi1
                real( trace( Vbar * (eye(N) - (vref*vref')) ) ) <= varpi2
        cvx_end

        % Diagnostics in bar-space
        MBW = full((Wbar+Wbar')/2); MBV = full((Vbar+Vbar')/2);
        eW  = eig(MBW); eV = eig(MBV);
        rW  = max(eW)/max(sum(eW), eps);
        rV  = max(eV)/max(sum(eV), eps);

        hist.obj(it)     = cvx_optval;
        hist.rankW(it)   = rW;
        hist.rankV(it)   = rV;
        hist.tW(it)      = full(tW);
        hist.tV(it)      = full(tV);
        hist.varpi1(it)  = full(varpi1);
        hist.varpi2(it)  = full(varpi2);
        hist.qos_lhs(it) = full(real(trace(En * Wbar)));

        if ~opts.verbose
            fprintf('[DC it %d] obj=%.3e | rankW=%.3f rankV=%.3f | tW=%.2e tV=%.2e | varpi=[%.2e %.2e] | QoS: %.2e>=%.2e\n',...
                it, hist.obj(it), rW, rV, hist.tW(it), hist.tV(it), hist.varpi1(it), hist.varpi2(it), hist.qos_lhs(it), rhs_qos_n);
        end

        % Update refs with principal eigenvectors for next linearization
        [wref,~] = eigs(MBW,1);
        [vref,~] = eigs(MBV,1);

        % Simple stopping: small relative obj change & small DC slacks
        if it >= 2
            rel_impr = abs(hist.obj(it) - obj_prev) / max(1, abs(obj_prev));
            if (rel_impr < opts.tol_rel) && (hist.varpi1(it)+hist.varpi2(it) < 1e-3)
                hist.obj = hist.obj(1:it);
                hist.rankW = hist.rankW(1:it);
                hist.rankV = hist.rankV(1:it);
                hist.tW = hist.tW(1:it);
                hist.tV = hist.tV(1:it);
                hist.varpi1 = hist.varpi1(1:it);
                hist.varpi2 = hist.varpi2(1:it);
                hist.qos_lhs = hist.qos_lhs(1:it);
                break;
            end
        end
        obj_prev = hist.obj(it);
    end

    % ------------------ unscale to original -----------------
    Wtilde = sW * Wbar;
    Vtilde = sV * Vbar;

    % Final principal eigenvectors (original scale is just scalar*MBW/MBV)
    [wmax,~] = eigs((Wtilde+Wtilde')/2, 1);
    [vmax,~] = eigs((Vtilde+Vtilde')/2, 1);

    % Feasibility diagnostics (in normalized/scaled coords)
    feas = struct();
    feas.status   = cvx_status;
    feas.qos_lhs  = hist.qos_lhs(end);
    feas.qos_rhs  = rhs_qos_n;
    feas.rankW    = hist.rankW(end);
    feas.rankV    = hist.rankV(end);
    feas.tW       = hist.tW(end);
    feas.tV       = hist.tV(end);
    feas.varpi1   = hist.varpi1(end);
    feas.varpi2   = hist.varpi2(end);
    feas.sW       = sW; feas.sV = sV; feas.tE = tE;
end