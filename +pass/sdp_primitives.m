function D = sdp_primitives(P, B, PT)
%SDP_PRIMITIVES  Build constants/matrices for the inner SDP (fixed positions).
%   Inputs:
%     P  - config (config_pass)
%     B  - beamformers struct from pass.beamformers(P,S): fields w, v, eta_vec
%     PT - transmit power (Watts)
%
%   Outputs (struct D):
%     E           - (N x N)  E = eta*eta'
%     Wx, Vx      - (N x N)  Wx = w*w', Vx = v*v'
%     gammaQoS    - scalar   (2^{R_QoS}-1)/PT
%     rmin_c/s    - scalars  lower-bound distances used for diag caps
%     diagW_cap   - (N x 1)  per-diagonal caps for Wtilde
%     diagV_cap   - (N x 1)  per-diagonal caps for Vtilde
%     wmax0/vmax0 - (N x 1)  principal eigenvectors of Wx, Vx (DC init)
%     N           - integer  number of elements

    N = numel(B.w);

    % Core matrices
    eta   = B.eta_vec;              % (N x 1)
    E     = eta * eta.';            % (N x N)
    Wx    = B.w * B.w';             % (N x N)
    Vx    = B.v * B.v';             % (N x N)

    % QoS parameter
    gammaQoS = (2^(P.R_QoS) - 1) / PT;

    % Lower-bound distances (paper’s definitions)
    rmin_c = sqrt( (P.rc*sin(P.phic))^2 + P.d^2 );
    rmin_s = sqrt( (P.rs*sin(P.phis))^2 + P.d^2 );

    % Diagonal caps (same value on all diagonals; keep as vectors)
    diagW_cap = (1/(N * rmin_c^2)) * ones(N,1);
    diagV_cap = (1/(N * rmin_s^2)) * ones(N,1);

    % Initial DC linearization points (top eigenvectors)
    [wmax0, ~] = eigs(Wx, 1);
    [vmax0, ~] = eigs(Vx, 1);

    % Pack
    D = struct( ...
        'E', E, 'Wx', Wx, 'Vx', Vx, ...
        'gammaQoS', gammaQoS, ...
        'rmin_c', rmin_c, 'rmin_s', rmin_s, ...
        'diagW_cap', diagW_cap, 'diagV_cap', diagV_cap, ...
        'wmax0', wmax0, 'vmax0', vmax0, ...
        'N', N );
end