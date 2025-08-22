function S = geom_channel(P, xp)
%GEOM_CHANNEL  Build geometry and channel primitives for PASS.
%   Inputs:
%     P  - config struct from config_pass()
%     xp - (N x 1) positions of pinching antennas along the x-axis (meters)
%
%   Outputs (struct S):
%     rc_n    - (N x 1) distances from each pinching element to the user
%     rs_n    - (N x 1) distances from each pinching element to the target
%     theta_n - (N x 1) in-waveguide radiation phase at each element
%     psi_p   - (N x 3) element coordinates [x, y, z]
%     psi0    - (1 x 3) feed-point coordinate
%     psi_c   - (1 x 3) user coordinate
%     psi_s   - (1 x 3) target coordinate

    N = numel(xp);
    yp = zeros(N,1);
    zp = P.d * ones(N,1);

    % Element & reference positions
    psi_p = [xp, yp, zp];                  % (N x 3)
    psi0  = [0, 0, P.d];                   % feed point (assumed on guide)
    psi_c = [P.rc*cos(P.phic), P.rc*sin(P.phic), 0];
    psi_s = [P.rs*cos(P.phis), P.rs*sin(P.phis), 0];

    % Near-field (spherical) distances to user/target
    rc_n = sqrt( sum( (psi_p - psi_c).^2, 2 ) );   % (N x 1)
    rs_n = sqrt( sum( (psi_p - psi_s).^2, 2 ) );   % (N x 1)

    % In-waveguide phase from feed to each element
    L_feed_n = sqrt( sum( (psi_p - psi0).^2, 2 ) );  % geometric distance
    theta_n  = 2*pi * P.eta_eff .* L_feed_n ./ P.lambda;

    % Pack outputs
    S = struct( ...
        'rc_n',    rc_n, ...
        'rs_n',    rs_n, ...
        'theta_n', theta_n, ...
        'psi_p',   psi_p, ...
        'psi0',    psi0, ...
        'psi_c',    psi_c, ...
        'psi_s',    psi_s );
end