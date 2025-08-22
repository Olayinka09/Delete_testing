function M = objective_metrics(P, B, PT)
%OBJECTIVE_METRICS  Compute R and Ps from w,v
%   Inputs:
%     P  - config struct (config_pass)
%     B  - struct from pass.beamformers(P,S): fields w, v, eta_vec
%     PT - transmit power (Watts)
%
%   Outputs (struct M):
%     R        - achievable rate (bps/Hz)
%     Ps       - illumination power at target (Watts, as modeled)
%     Ps_dBW   - 10*log10(Ps)
%     Ps_dBm   - 10*log10(Ps) + 30
%     qos_ok   - logical: R >= P.R_QoS

    eta = B.eta_vec;
    w = B.w; v = B.v;

    M.R      = log2( 1 + PT * abs(eta' * w)^2 / P.sigma2 );
    M.Ps     = real( PT * (eta' * v) * (v' * eta) );  % scalar
    M.Ps_dBW = 10*log10(M.Ps);
    M.Ps_dBm = M.Ps_dBW + 30;
    M.qos_ok = (M.R >= P.R_QoS);
end