function [Wtilde, Vtilde, obj, feas] = inner_sdp_penalty_only(P, D, rho)
%INNER_SDP_PENALTY_ONLY  Convex inner SDP with Frobenius penalties.
% Fixed positions x^p. No DC (rank) surrogate yet.
%
% Inputs:
%   P   : config_pass()
%   D   : pass.sdp_primitives(P,B,PT) (E, Wx, Vx, gammaQoS, diag caps, N)
%   rho : penalty weight (e.g., 1 or 10)
%
% Outputs:
%   Wtilde, Vtilde : (N x N) Hermitian PSD in ORIGINAL scale
%   obj            : objective value
%   feas           : diagnostics

    N = D.N;

    % ---- scale diagonal caps to ~1
    sW = D.diagW_cap(1);    % >0
    sV = D.diagV_cap(1);
    if sW <= 0 || sV <= 0, error('Diagonal caps must be positive.'); end

    % bar-space variables: Wbar = Wtilde/sW, Vbar = Vtilde/sV
    rhs_qos  = D.gammaQoS * P.sigma2 / sW;

    % normalize E by its trace to reduce condition number
    tE  = trace(D.E) + eps;
    En  = D.E / tE;
    rhs_qos_n = rhs_qos / tE;

    % targets in bar-space for penalties
    Wbar0 = D.Wx / sW;   % complex Hermitian
    Vbar0 = D.Vx / sV;

    cvx_clear
    try
        cvx_solver sedumi
    catch
        cvx_solver sdpt3
        cvx_solver_settings('rmdepconstr', 1)
    end
    cvx_precision best

    cvx_begin sdp
        % Hermitian (complex) PSD variables
        variable Wbar(N,N) hermitian semidefinite
        variable Vbar(N,N) hermitian semidefinite
        % Epigraph variables for Frobenius norms
        variable tW nonnegative
        variable tV nonnegative

        % Objective: sensing term - Frobenius penalties
        maximize( sV * trace(En * Vbar) - (1/(2*rho))*(tW + tV) )

        subject to
            % Epigraph constraints
            norm(Wbar - Wbar0, 'fro') <= tW
            norm(Vbar - Vbar0, 'fro') <= tV

            % Diagonal caps (scaled)
            diag(Wbar) <= 1
            diag(Vbar) <= 1

            % QoS (scaled & normalized)
            trace(En * Wbar) >= rhs_qos_n
    cvx_end

    % Unscale back to original space
    Wtilde = sW * Wbar;
    Vtilde = sV * Vbar;

    obj = cvx_optval;

    feas = struct();
    feas.status   = cvx_status;
    feas.qos_lhs  = full(real(trace(En * Wbar)));
    feas.qos_rhs  = rhs_qos_n;
    feas.tW       = full(tW);
    feas.tV       = full(tV);
    feas.diagW    = full(diag(Wbar));
    feas.diagV    = full(diag(Vbar));
    feas.sW       = sW;
    feas.sV       = sV;
    feas.tE       = tE;
end
