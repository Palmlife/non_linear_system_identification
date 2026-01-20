clear; close all; clc;
fs = 4000;
fmax = 1200;  % Max excited frequency

%% Non-parametric estimation
[G_ML, f, noise_var, total_var] = robust_method_SA("results_2/fullMOut", fs, 1);

%% Parametric estimation
na = 5;
nb = 5;

%% Parametric estimation with optimal order - Compare methods
% Use total variance (noise + stochastic distortion) for weighted methods
[b_LLS, a_LLS, G_LLS]  = parametric_LS(G_ML, f, fs, nb, na, [], 'LLS', fmax);
[b_TLS, a_TLS, G_TLS]  = parametric_LS(G_ML, f, fs, nb, na, [], 'TLS', fmax);
[b_GTLS, a_GTLS, G_GTLS] = parametric_LS(G_ML, f, fs, nb, na, total_var, 'GTLS', fmax);
[b_WLS, a_WLS, G_WLS]  = parametric_LS(G_ML, f, fs, nb, na, total_var, 'WLS', fmax);
[b_IWLS, a_IWLS, G_IWLS, n_iter_IWLS] = parametric_LS(G_ML, f, fs, nb, na, [], 'IWLS', fmax);
[b_IQML, a_IQML, G_IQML, n_iter_IQML] = parametric_LS(G_ML, f, fs, nb, na, total_var, 'IQML', fmax);
[b_ML, a_ML, G_ML_est, n_iter_ML] = parametric_LS(G_ML, f, fs, nb, na, total_var, 'ML', fmax);

% Print transfer functions for LaTeX
fprintf('\n=== Estimated Transfer Functions (na=%d, nb=%d) ===\n', na, nb);
print_tf('LLS', b_LLS, a_LLS);
print_tf('TLS', b_TLS, a_TLS);
print_tf('GTLS', b_GTLS, a_GTLS);
print_tf('WLS', b_WLS, a_WLS);
print_tf('IWLS', b_IWLS, a_IWLS);
print_tf('IQML', b_IQML, a_IQML);
print_tf('ML', b_ML, a_ML);

% Compute weighted and unweighted costs
idx = (f > 0) & (f < fmax);
cost_weighted = @(G) sum(abs(G_ML(idx) - G(idx)).^2 ./ total_var(idx));
cost_unweighted = @(G) sum(abs(G_ML(idx) - G(idx)).^2);
fprintf('\nWeighted costs - LLS: %.2e, TLS: %.2e, GTLS: %.2e, WLS: %.2e, IWLS: %.2e, IQML: %.2e, ML: %.2e\n', ...
    cost_weighted(G_LLS), cost_weighted(G_TLS), cost_weighted(G_GTLS), cost_weighted(G_WLS), cost_weighted(G_IWLS), cost_weighted(G_IQML), cost_weighted(G_ML_est));
fprintf('Unweighted costs - LLS: %.2e, TLS: %.2e, GTLS: %.2e, WLS: %.2e, IWLS: %.2e, IQML: %.2e, ML: %.2e\n', ...
    cost_unweighted(G_LLS), cost_unweighted(G_TLS), cost_unweighted(G_GTLS), cost_unweighted(G_WLS), cost_unweighted(G_IWLS), cost_unweighted(G_IQML), cost_unweighted(G_ML_est));
fprintf('Iterations - IWLS: %d, IQML: %d, ML: %d\n', n_iter_IWLS, n_iter_IQML, n_iter_ML);

%% Generate figures for report
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultLineLineWidth', 1);

% Use absolute path (works in VS Code MATLAB)
scriptDir = fileparts(mfilename('fullpath'));
picDir = fullfile(scriptDir, '..', 'report', 'pic');

methods = {'LLS', 'TLS', 'GTLS', 'WLS', sprintf('IWLS (n=%d)', n_iter_IWLS), sprintf('IQML (n=%d)', n_iter_IQML), sprintf('ML (n=%d)', n_iter_ML)};
methods_short = {'LLS', 'TLS', 'GTLS', 'WLS', 'IWLS', 'IQML', 'ML'};  % For filenames
G_est = {G_LLS, G_TLS, G_GTLS, G_WLS, G_IWLS, G_IQML, G_ML_est};
colors = {'b', 'r', 'c', 'm', [0.2 0.6 0.2], 'g', [0.8 0.4 0]};

for i = 1:7
    fig = figure('Units', 'centimeters', 'Position', [5 5 14 10]);
    
    % Magnitude
    subplot(2,1,1);
    plot(f, db(G_ML), 'k.', 'MarkerSize', 2); hold on;
    plot(f, db(G_est{i}), '-', 'Color', colors{i}, 'LineWidth', 1);
    ylabel('Magnitude [dB]'); xlim([0 fmax]); grid on;
    legend('Measured', methods{i}, 'Location', 'best');
    
    % Residuals
    subplot(2,1,2);
    plot(f(idx), db(abs(G_ML(idx)-G_est{i}(idx))), '.', 'Color', colors{i}, 'MarkerSize', 2); hold on;
    plot(f(idx), db(sqrt(total_var(idx))), 'k-', 'LineWidth', 1);
    xlabel('Frequency [Hz]'); ylabel('[dB]'); xlim([0 fmax]); grid on;
    %the total standard deviation is plotted as a reference
    legend('Residual', 'Total std', 'Location', 'best');
    
    exportgraphics(fig, fullfile(picDir, sprintf('2_%s_Fit.png', methods_short{i})), 'Resolution', 300);
end

% Weighted cost comparison bar chart
fig = figure('Units', 'centimeters', 'Position', [5 5 14 7]);
bar([cost_weighted(G_LLS), cost_weighted(G_TLS), cost_weighted(G_GTLS), cost_weighted(G_WLS), cost_weighted(G_IWLS), cost_weighted(G_IQML), cost_weighted(G_ML_est)], 'FaceColor', [0.3 0.5 0.7]);
set(gca, 'XTickLabel', methods_short);
ylabel('Weighted cost'); grid on;
exportgraphics(fig, fullfile(picDir, '2_Cost_Comparison.png'), 'Resolution', 300);

% Unweighted cost comparison bar chart
fig = figure('Units', 'centimeters', 'Position', [5 5 14 7]);
bar([cost_unweighted(G_LLS), cost_unweighted(G_TLS), cost_unweighted(G_GTLS), cost_unweighted(G_WLS), cost_unweighted(G_IWLS), cost_unweighted(G_IQML), cost_unweighted(G_ML_est)], 'FaceColor', [0.7 0.3 0.3]);
set(gca, 'XTickLabel', methods_short);
ylabel('Unweighted cost'); grid on;
exportgraphics(fig, fullfile(picDir, '2_Cost_Comparison_Unweighted.png'), 'Resolution', 300);

% Weighted cost comparison bar chart WITHOUT IWLS (outlier)
fig = figure('Units', 'centimeters', 'Position', [5 5 14 7]);
methods_no_iwls = {'LLS', 'TLS', 'GTLS', 'WLS', 'IQML', 'ML'};
bar([cost_weighted(G_LLS), cost_weighted(G_TLS), cost_weighted(G_GTLS), cost_weighted(G_WLS), cost_weighted(G_IQML), cost_weighted(G_ML_est)], 'FaceColor', [0.3 0.5 0.7]);
set(gca, 'XTickLabel', methods_no_iwls);
ylabel('Weighted cost'); grid on;
exportgraphics(fig, fullfile(picDir, '2_Cost_Comparison_NoIWLS.png'), 'Resolution', 300);

% Unweighted cost comparison bar chart WITHOUT WLS (outlier)
fig = figure('Units', 'centimeters', 'Position', [5 5 14 7]);
methods_no_wls = {'LLS', 'TLS', 'GTLS', 'IWLS', 'IQML', 'ML'};
bar([cost_unweighted(G_LLS), cost_unweighted(G_TLS), cost_unweighted(G_GTLS), cost_unweighted(G_IWLS), cost_unweighted(G_IQML), cost_unweighted(G_ML_est)], 'FaceColor', [0.7 0.3 0.3]);
set(gca, 'XTickLabel', methods_no_wls);
ylabel('Unweighted cost'); grid on;
exportgraphics(fig, fullfile(picDir, '2_Cost_Comparison_Unweighted_NoWLS.png'), 'Resolution', 300);

fprintf('Figures saved to report/pic/\n');

%% Helper function to print transfer function
function print_tf(name, b, a)
    fprintf('\n%s:\n', name);
    fprintf('B(z) = ');
    for i = 1:length(b)
        if i == 1
            fprintf('%.4e', b(i));
        else
            if b(i) >= 0
                fprintf(' + %.4e z^{-%d}', b(i), i-1);
            else
                fprintf(' - %.4e z^{-%d}', abs(b(i)), i-1);
            end
        end
    end
    fprintf('\n');
    fprintf('A(z) = 1');
    for i = 2:length(a)
        if a(i) >= 0
            fprintf(' + %.4e z^{-%d}', a(i), i-1);
        else
            fprintf(' - %.4e z^{-%d}', abs(a(i)), i-1);
        end
    end
    fprintf('\n');
end