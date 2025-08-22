function B = beamformers(P, S)
%BEAMFORMERS  Build pinching beamformers w(x^p) and v(x^p)
%   Inputs:
%     P - config struct (config_pass)
%     S - struct from pass.geom_channel(P, xp), must contain:
%         rc_n, rs_n, theta_n
%
%   Outputs:
%     B.w       - (N x 1) communication beamformer w(x^p)
%     B.v       - (N x 1) sensing beamformer v(x^p)
%     B.alpha_n - (N x 1) per-element power coefficients (equal power)
%     B.eta_vec - (N x 1) constant vector eta^{1/2}

    N = numel(S.theta_n);

    % Equal power allocation (alpha_s split evenly)
    alpha_n = (P.alpha_sum / P.N) * ones(N,1);

    % Phases toward user/target including in-waveguide phase
    phase_c = (2*pi/P.lambda) .* S.rc_n + S.theta_n;
    phase_s = (2*pi/P.lambda) .* S.rs_n + S.theta_n;

    % Beamformers (eqs. (5)–(6))
    w = exp(-1j .* phase_c) ./ (sqrt(alpha_n) .* S.rc_n);
    v = exp(-1j .* phase_s) ./ (sqrt(alpha_n) .* S.rs_n);

    % Constant vector eta^{1/2}
    eta_vec = sqrt(P.eta) * ones(N,1);

    % Pack
    B = struct('w', w, 'v', v, 'alpha_n', alpha_n, 'eta_vec', eta_vec);
end