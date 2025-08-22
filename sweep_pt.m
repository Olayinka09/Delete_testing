function T = sweep_pt()
clc;

% --- solver set once in console:
% try, cvx_end; end; cvx_clear; cvx_solver sedumi; cvx_precision best

P    = config_pass();
xp0  = linspace(-P.L/2, P.L/2, P.N).';
PTdBm_list = 65:5:90;

% Alg.1 / Alg.2 settings (paper-faithful)
optsA1 = struct('rho',10,'rho_i',1,'c1bar',0.8,'eps1',1e-3,'max_dc_iter',20,'verbose',false);
rho0=10; c2=0.8; eps2=1e-3; eps3=1e-2;
optsA2 = struct('dc_iter',3,'rho_i',1,'c1bar',0.8,'eps1',1e-3,...
                'max_outer',10,'max_inner',20,'verbose',true);

K = numel(PTdBm_list);
R_base  = zeros(K,1);  Ps_base_dBm = zeros(K,1);  QoS_b = false(K,1);
R_opt   = zeros(K,1);  Ps_opt_dBm  = zeros(K,1);  QoS_o = false(K,1);

for k = 1:K
    PTdBm = PTdBm_list(k);
    PT    = 10^((PTdBm-30)/10);

    % ---- Baseline @ uniform layout (Alg.1 once)
    S0 = pass.geom_channel(P, xp0);
    B0 = pass.beamformers(P, S0);
    D0 = pass.sdp_primitives(P, B0, PT);
    [~, ~, wb, vb, ~, ~] = pass.inner_sdp_dc_solve(P, D0, optsA1);
    phi0 = -angle(wb(1)); wb = wb*exp(1j*phi0); vb = vb*exp(1j*phi0);
    B0.w = wb; B0.v = vb;
    M0 = pass.objective_metrics(P, B0, PT);

    % ---- Optimized PASS (Alg.2) starting from same xp0
    [xp_opt, ~] = pass.ao_run_paper(P, xp0, PT, rho0, c2, eps2, eps3, optsA2);
    S1 = pass.geom_channel(P, xp_opt);
    B1 = pass.beamformers(P, S1);
    D1 = pass.sdp_primitives(P, B1, PT);
    [~, ~, wo, vo, ~, ~] = pass.inner_sdp_dc_solve(P, D1, optsA1);
    phi1 = -angle(wo(1)); wo = wo*exp(1j*phi1); vo = vo*exp(1j*phi1);
    B1.w = wo; B1.v = vo;
    M1 = pass.objective_metrics(P, B1, PT);

    % store
    R_base(k)     = M0.R;   Ps_base_dBm(k) = M0.Ps_dBm;  QoS_b(k) = M0.qos_ok;
    R_opt(k)      = M1.R;   Ps_opt_dBm(k)  = M1.Ps_dBm;  QoS_o(k) = M1.qos_ok;

    fprintf('[%ddBm] Base: R=%.3f | Ps=%.3f dBm | QoS=%d  ||  Opt: R=%.3f | Ps=%.3f dBm | QoS=%d\n',...
        PTdBm, M0.R, M0.Ps_dBm, M0.qos_ok, M1.R, M1.Ps_dBm, M1.qos_ok);
end

% Table of results
T = table(PTdBm_list(:), R_base, R_opt, Ps_base_dBm, Ps_opt_dBm, QoS_b, QoS_o, ...
          'VariableNames', {'PT_dBm','R_base','R_opt','Ps_base_dBm','Ps_opt_dBm','QoS_base','QoS_opt'});

% Quick plot like Fig. 2 (illumination vs PT)
figure; hold on; grid on;
plot(PTdBm_list, Ps_opt_dBm,  'o-','LineWidth',1.5);
plot(PTdBm_list, Ps_base_dBm, 's-','LineWidth',1.5);
xlabel('Transmit power at BS (dBm)');
ylabel('Illumination power (dBm)');
legend({'Pinching antenna (optimized)','Fixed pinching (uniform)'},'Location','NorthWest');
title('Illumination power vs. transmit power');

end