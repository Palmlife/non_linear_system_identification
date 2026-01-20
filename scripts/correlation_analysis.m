% Correlation analysis between |A(z)|^2 and total variance
% This script analyzes why LLS performs well despite not using noise information
clear; close all; clc;

fs = 4000;
fmax = 1200;  % Max excited frequency

%% Non-parametric estimation
[G_ML, f, noise_var, total_var] = robust_method_SA("results_2/fullMOut", fs, 1);

%% Parametric estimation with LLS
na = 5;
nb = 5;
[b_LLS, a_LLS, G_LLS] = parametric_LS(G_ML, f, fs, nb, na, [], 'LLS', fmax);

%% Analyze correlation between |A|^2 and total_var
idx = (f > 0) & (f < fmax);
z = exp(1j*2*pi*f/fs);
A_LLS = polyval(flip(a_LLS), z.^-1);

% Compute correlation
corr_coef = corrcoef(log(total_var(idx)), log(abs(A_LLS(idx)).^2));
fprintf('Correlation between log(total_var) and log(|A_LLS|^2): %.4f\n', corr_coef(1,2));

% Variance statistics
fprintf('\nTotal variance stats:\n');
fprintf('  min: %.2e, max: %.2e, ratio: %.1f\n', ...
    min(total_var(idx)), max(total_var(idx)), max(total_var(idx))/min(total_var(idx)));

fprintf('\n|A|^2 stats:\n');
fprintf('  min: %.2e, max: %.2e, ratio: %.1f\n', ...
    min(abs(A_LLS(idx)).^2), max(abs(A_LLS(idx)).^2), max(abs(A_LLS(idx)).^2)/min(abs(A_LLS(idx)).^2));

%% Generate figure
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultLineLineWidth', 1);

scriptDir = fileparts(mfilename('fullpath'));
picDir = fullfile(scriptDir, '..', 'report', 'pic');

fig = figure('Units', 'centimeters', 'Position', [5 5 14 10]);

% Subplot 1: Both quantities vs frequency
subplot(2,1,1);
yyaxis left
semilogy(f(idx), abs(A_LLS(idx)).^2, 'b-', 'LineWidth', 1);
ylabel('|A(z)|^2');
yyaxis right
semilogy(f(idx), total_var(idx), 'r-', 'LineWidth', 1);
ylabel('\sigma^2 (total variance)');
xlabel('Frequency [Hz]');
title(sprintf('Correlation: %.2f', corr_coef(1,2)));
grid on;
legend('|A(z)|^2', '\sigma^2', 'Location', 'best');

% Subplot 2: Scatter plot
subplot(2,1,2);
loglog(abs(A_LLS(idx)).^2, total_var(idx), 'k.', 'MarkerSize', 3);
hold on;
% Fit line
p = polyfit(log(abs(A_LLS(idx)).^2), log(total_var(idx)), 1);
x_fit = logspace(log10(min(abs(A_LLS(idx)).^2)), log10(max(abs(A_LLS(idx)).^2)), 100);
y_fit = exp(polyval(p, log(x_fit)));
loglog(x_fit, y_fit, 'r-', 'LineWidth', 1.5);
xlabel('|A(z)|^2'); ylabel('\sigma^2 (total variance)');
title(sprintf('Log-log slope: %.2f (corr = %.2f)', p(1), corr_coef(1,2)));
legend('Data', sprintf('Fit: \\sigma^2 \\propto |A|^{%.1f}', 2*p(1)), 'Location', 'best');
grid on;

exportgraphics(fig, fullfile(picDir, '2_A_vs_Variance_Correlation.png'), 'Resolution', 300);
fprintf('\nFigure saved to report/pic/2_A_vs_Variance_Correlation.png\n');
