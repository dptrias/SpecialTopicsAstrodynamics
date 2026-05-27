%% WP2
script_dir = fileparts(mfilename('fullpath'));
% System parameters
mu_SE = 3.0542e-6;
mu_SM = 3.227154996101724e-7;

% Calculate Lagrange points for Sun-Earth and Sun-Mars systems
L2_SE = ClassicalLagrangePoints(mu_SE, 2);
L2_SM = ClassicalLagrangePoints(mu_SM, 2);

dx = 0.00050;  % Desired shift in the Lagrange point (in LU)
beta_SE = beta_ArtificialEqPoint(L2_SE - dx, mu_SE);
beta_SM = beta_ArtificialEqPoint(L2_SM - dx, mu_SM);

fprintf('Beta for Sun-Earth L2: %.8e\n', beta_SE);
fprintf('Beta for Sun-Mars L2: %.8e\n', beta_SM);

mu_S = 132712e15;
m_s = 100;
S_Sun = 1361;  % W/m^2
c = 299792458;  % m/s
LU_SE = 149597870.7e3;  % m
LU_SM = 208321282e3;  % m
P_Earth = S_Sun / c;  % Radiation pressure at Earth distance
P_Mars = S_Sun * (LU_SE / LU_SM)^2 / c;  % Radiation pressure at Mars distance

x_SE = (L2_SE - dx) * LU_SE;
x_SM = (L2_SM - dx) * LU_SM;
A_SE = beta_SE * mu_S * m_s / (2 * SolarRadiationPressure(S_Sun, x_SE) * x_SE^2);
A_SM = beta_SM * mu_S * m_s / (2 * SolarRadiationPressure(S_Sun, x_SM) * x_SM^2);

fprintf('Required length for Sun-Earth L2 sail: %.4f m\n', sqrt(A_SE));
fprintf('Required length for Sun-Mars L2 sail: %.4f m\n', sqrt(A_SM));

% Sensitivity analysis for shift in the Lagrange point
dx_values = linspace(dx, 2*dx, 21);
beta_SE_values = arrayfun(@(dx_i) beta_ArtificialEqPoint(L2_SE - dx_i, mu_SE), dx_values);
beta_SM_values = arrayfun(@(dx_i) beta_ArtificialEqPoint(L2_SM - dx_i, mu_SM), dx_values);

figure_style;
figure;
plot(dx_values, beta_SE_values, 'b-', 'LineWidth', 2, 'DisplayName', 'Sun-Earth L2');
hold on;
plot(dx_values, beta_SM_values, 'r-', 'LineWidth', 2, 'DisplayName', 'Sun-Mars L2');
xlabel('Shift in Lagrange point, $\Delta x$ [LU]');
ylabel('Lightness number, $\beta$ [-]');
%title('Sensitivity of Lightness Number to Shift in Lagrange Point');
xlim([dx, 2*dx]);
legend('Location', 'northwest');
grid on;
saveas(gcf, fullfile(script_dir, '../figures', 'beta_sensitivity.png'));

beta_test = beta_ArtificialEqPoint(0.983867, mu_SE);
fprintf('Calculated beta for Sun-Earth AEP at x=0.983867 LU: %.8e\n', beta_test);

function beta = beta_ArtificialEqPoint(x_AEP, mu)
    beta = (x_AEP + mu)^2 * (- x_AEP + (1 - mu) / (x_AEP + mu)^2 + mu / (x_AEP - (1 - mu))^2) / (1 - mu);
end

function p = SolarRadiationPressure(S, r)
    c = 299792458;  % Speed of light in m/s
    p = S / c * (149597870.7e3 / r)^2;  % Radiation pressure at distance r from the Sun
end