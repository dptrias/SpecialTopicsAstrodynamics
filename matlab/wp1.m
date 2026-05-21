%% WP1
% Earth-Moon system parameters
mu_SE = 3.0542e-6;
LU_SE = 149597870.7; % km
TU_SE = 365.256 * 24 * 3600; % s
% Sun-Mars system parameters
mu_SM = 3.227154996101724e-7;
LU_SM = 208321282; % km
TU_SM = 8253622; % s

% Calculate Lagrange points for Sun-Earth and Sun-Mars systems
L2_SE = ClassicalLagrangePoints(mu_SE, 2);
L2_SM = ClassicalLagrangePoints(mu_SM, 2);

L2_SE_JPL = 1.01009044;
L2_SM_JPL = 1.00476311;

fprintf('Sun-Earth L2: %.8f (JPL: %.8f)\n', L2_SE, L2_SE_JPL);
fprintf('Sun-Mars L2: %.8f (JPL: %.8f)\n', L2_SM, L2_SM_JPL);

fprintf('Difference for Sun-Earth L2: %.8e\n', abs(L2_SE - L2_SE_JPL));
fprintf('Difference for Sun-Mars L2: %.8e\n', abs(L2_SM - L2_SM_JPL));