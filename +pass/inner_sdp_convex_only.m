function [Wtilde, Vtilde, obj, feas] = inner_sdp_convex_only(P, D)
%INNER_SDP_CONVEX_ONLY  Convex-only inner SDP for fixed positions x^p.
%   No DC (rank) surrogate and no Frobenius penalties — just the core
%   PSD, diagonal caps, and QoS constraint. Numerically stabilized via
%   variable scaling and normalization. Results are unscaled on output.
%
% Inputs:
%   P  : config struct from config_pass()
%   D  : struct from pass.sdp_primitives(P,B,PT) with fields:
%        E, gammaQoS, diagW_cap, diagV_cap, N, and P.sigma2 via P
%
% Outputs:
%   Wtilde, Vtilde : (N x N) Hermitian PSD matrices in ORIGINAL scale
%   obj            : CVX objective value (scalar)
%   feas           : diagnostics (status, QoS lhs/rhs in normalized scale, etc.)

    N = D.N;

    % --- Scale so diagonal caps become ~1 -------------------------------
    % We define Wbar = Wtilde / sW, Vbar = Vtilde / sV so that:
    %   diag(Wbar) <= 1, diag(Vbar) <= 1.
    % (diag caps are equal across n by construction.)
    sW = D.diagW_cap(1);    % > 0
    sV = D.diagV_cap(1);    % > 0
    if sW <= 0 || sV <= 0
        error('Diagonal caps must be positive.');
    end

    % QoS right-hand side in the scaled (bar) variables:
    % trace(E * (sW*Wbar)) >= gammaQoS*sigma2  =>  trace(E*Wbar) >= (gammaQoS*sigma2)/sW
    rhs_qos = D.gammaQoS * P.sigma2 / sW;

    % --- Normalize the QoS matrix E to reduce condition number ----------
    % Let En = E / trace(E). Apply the same normalization to rhs.
    tE  = trace(D.E) + eps;
    En  = D.E / tE;
    rhs_qos_n = rhs_qos / tE;

    % --- Solve the stabilized convex SDP --------------------------------
    cvx_clear
    % Prefer SeDuMi for this stabilized form; fall back to SDPT3 with dep-removal
    try
        cvx_solver sedumi
    catch
        cvx_solver sdpt3
        cvx_solver_settings('rmdepconstr', 1)
    end
    cvx_precision best

    cvx_begin sdp
        variable Wbar(N,N) symmetric semidefinite
        variable Vbar(N,N) symmetric semidefinite

        % Objective: maximize sensing term (scaled). A positive scalar sV
        % simply scales the objective; it does not affect the optimizer.
        maximize( sV * trace(En * Vbar) )

        subject to
            % Diagonal caps (now ~1)
            diag(Wbar) <= 1;
            diag(Vbar) <= 1;

            % QoS constraint in normalized/scaled coordinates
            trace(En * Wbar) >= rhs_qos_n;
    cvx_end

    % --- Unscale back to original variables -----------------------------
    Wtilde = sW * Wbar;
    Vtilde = sV * Vbar;

    obj = cvx_optval;

    % Diagnostics (note: QoS reported in normalized/scaled coordinates)
    feas = struct();
    feas.status   = cvx_status;
    feas.qos_lhs  = full(real(trace(En * Wbar)));  % compare with rhs_qos_n
    feas.qos_rhs  = rhs_qos_n;
    feas.diagW    = full(diag(Wbar));
    feas.diagV    = full(diag(Vbar));
    feas.sW       = sW;
    feas.sV       = sV;
    feas.tE       = tE;
end