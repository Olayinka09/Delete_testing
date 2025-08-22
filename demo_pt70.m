clc;
try, cvx_end; end
cvx_clear

% 1) Config & uniform layout
P    = config_pass();
xp0  = linspace(-P.L/2, P.L/2, P.N).';
PTdBm = 65; 
PT    = 10^((PTdBm-30)/10);

% 2) BASELINE @ xp0: Alg.1 -> (w,v) -> metrics
S0 = pass.geom_channel(P, xp0);
B0 = pass.beamformers(P, S0);
D0 = pass.sdp_primitives(P, B0, PT);

optsA1 = struct('rho',10,'rho_i',1,'c1bar',0.8,'eps1',1e-3,'max_dc_iter',20,'verbose',false);
[Wb, Vb, wb, vb, ~, feas0] = pass.inner_sdp_dc_solve(P, D0, optsA1);

% phase-align beams and put into B0
phi0 = -angle(wb(1)); 
wb   = wb * exp(1j*phi0); 
vb   = vb * exp(1j*phi0);
B0.w = wb; 
B0.v = vb;

M_base = pass.objective_metrics(P, B0, PT);

% 3) OPTIMIZED (Algorithm 2), then Alg.1 once at final layout -> metrics
rho0=10; c2=0.8; eps2=1e-3; eps3=1e-2;  % bar-space gap for Alg.2
optsA2 = struct('dc_iter',3,'rho_i',1,'c1bar',0.8,'eps1',1e-3,...
                'max_outer',10,'max_inner',20,'verbose',true);

[xp_opt, out] = pass.ao_run_paper(P, xp0, PT, rho0, c2, eps2, eps3, optsA2);

S1 = pass.geom_channel(P, xp_opt);
B1 = pass.beamformers(P, S1);
D1 = pass.sdp_primitives(P, B1, PT);

[Wo, Vo, wo, vo, ~, feas1] = pass.inner_sdp_dc_solve(P, D1, optsA1);

phi1 = -angle(wo(1)); 
wo   = wo * exp(1j*phi1); 
vo   = vo * exp(1j*phi1);
B1.w = wo; 
B1.v = vo;

M_opt = pass.objective_metrics(P, B1, PT);

% 4) Print
fprintf('\n=== Metrics @ PT = %d dBm ===\n', PTdBm);
fprintf('Baseline (uniform): R = %.3f bps/Hz | Ps = %.3f dBm (%.3f dBW) | QoS=%d\n', ...
        M_base.R, M_base.Ps_dBm, M_base.Ps_dBW, M_base.qos_ok);
fprintf('Optimized (Alg.2):  R = %.3f bps/Hz | Ps = %.3f dBm (%.3f dBW) | QoS=%d\n', ...
        M_opt.R,  M_opt.Ps_dBm,  M_opt.Ps_dBW,  M_opt.qos_ok);