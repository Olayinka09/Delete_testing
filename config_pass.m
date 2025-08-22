function P = config_pass()
    % Defaults follow the paper's Numerical Results section.
    c = physconst('LightSpeed');

    % Core parameters
    P.fc     = 28e9;                 % carrier (Hz)
    P.lambda = c / P.fc;
    P.eta    = c^2 / (16*pi^2*P.fc^2); % path-loss constant used in paper vectors
    P.sigma2_dBm = -105;            % noise power
    P.sigma2 = 10^((P.sigma2_dBm - 30)/10);  % Watts

    % Geometry (waveguides along x-axis, height d)
    P.d   = 10;                     % height (m)
    P.L   = 50;                     % waveguide length (m)
    P.N   = 16;                     % number of pinching antennas (TX waveguide)
    P.DELTAx = 0.5;                 % min spacing (m) (placeholder for later)

    % User & target (x‑y plane)
    P.rs   = 30;                    % target range (m)
    P.phis = pi/3;                  % target azimuth (rad)
    P.rc   = 15*sqrt(2);            % user range (m)
    P.phic = 5*pi/4;                % user azimuth (rad)

    % In‑waveguide effective index
    P.eta_eff = 1.4;

    % QoS
    P.R_QoS = 10;                   % bps/Hz

    % Power allocation (equal power model)
    P.alpha_sum = 1;                % total radiation coefficient (alpha_s)
end