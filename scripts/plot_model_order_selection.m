% Model order selection comparison script
% Note: Iterative methods (IWLS, IQML, ML) excluded - too slow for grid search
clear; close all; clc;

fs = 4000;
fmax = 1200;
na_max = 10; nb_max = 10;

% Load non-parametric estimate
[G_ML, f, noise_var] = robust_method_SA("results_2/fullMOut", fs, 1);

% Compute for non-iterative methods only (fast)
methods = {'LLS', 'TLS', 'GTLS', 'WLS'};
n_methods = length(methods);
AIC_all = nan(na_max, nb_max+1, n_methods);
GIC_all = nan(na_max, nb_max+1, n_methods);
whiteness_all = nan(na_max, nb_max+1, n_methods);

parfor (m = 1:n_methods, n_methods)
    fprintf('Computing model order selection for %s...\n', methods{m});
    [~, ~, AIC_all(:,:,m), whiteness_all(:,:,m), GIC_all(:,:,m)] = model_order_selection(G_ML, f, fs, noise_var, na_max, nb_max, fmax, methods{m});
end

% Output directory
scriptDir = fileparts(mfilename('fullpath'));
picDir = fullfile(scriptDir, '..', 'report', 'pic');

% Create individual AICc + GIC + Whiteness plot for each method (3 subplots)
for m = 1:n_methods
    AIC = AIC_all(:,:,m);
    GIC = GIC_all(:,:,m);
    whiteness = whiteness_all(:,:,m);
    
    % Find AICc-optimal order
    [~, i_aic] = min(AIC(:)); 
    [na_aic, nb_aic] = ind2sub(size(AIC), i_aic); 
    nb_aic_plot = nb_aic - 1;
    
    % Find GIC-optimal order
    [~, i_gic] = min(GIC(:));
    [na_gic, nb_gic] = ind2sub(size(GIC), i_gic);
    nb_gic_plot = nb_gic - 1;
    
    fig = figure('Units', 'centimeters', 'Position', [5 5 20 6]);
    
    % AICc subplot
    subplot(131); 
    imagesc(0:nb_max, 1:na_max, AIC); set(gca,'YDir','normal');
    colorbar; xlabel('n_b'); ylabel('n_a'); 
    title(sprintf('AICc (%s)', methods{m})); hold on;
    plot(nb_aic_plot, na_aic, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    
    % GIC subplot
    subplot(132); 
    imagesc(0:nb_max, 1:na_max, GIC); set(gca,'YDir','normal');
    colorbar; xlabel('n_b'); ylabel('n_a'); 
    title(sprintf('GIC (%s)', methods{m})); hold on;
    plot(nb_gic_plot, na_gic, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Whiteness subplot
    subplot(133); 
    imagesc(0:nb_max, 1:na_max, whiteness*100); set(gca,'YDir','normal');
    colorbar; xlabel('n_b'); ylabel('n_a'); 
    title(sprintf('Whiteness %% (%s)', methods{m})); hold on;
    % Mark both AICc and GIC optimal points
    plot(nb_aic_plot, na_aic, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(nb_gic_plot, na_gic, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
    
    exportgraphics(fig, fullfile(picDir, sprintf('3_Model_Order_%s.png', methods{m})), 'Resolution', 300);
    close(fig);
end

% Summary comparison plot: All AICc in one figure (4x2 grid)
fig2 = figure('Units', 'centimeters', 'Position', [5 5 20 18]);
for m = 1:n_methods
    subplot(4,2,m); 
    imagesc(0:nb_max, 1:na_max, AIC_all(:,:,m)); set(gca,'YDir','normal');
    colorbar; xlabel('n_b'); ylabel('n_a'); 
    title(sprintf('AICc (%s)', methods{m})); hold on;
    [~, im] = min(AIC_all(:,:,m), [], 'all', 'linear');
    [na_m, nb_m] = ind2sub([na_max, nb_max+1], im);
    plot(nb_m-1, na_m, 'go', 'MarkerSize', 8, 'LineWidth', 2);
end
exportgraphics(fig2, fullfile(picDir, '3_AICc_All_Methods.png'), 'Resolution', 300);

% Summary comparison plot: All GIC in one figure (4x2 grid)
fig2b = figure('Units', 'centimeters', 'Position', [5 5 20 18]);
for m = 1:n_methods
    subplot(4,2,m); 
    imagesc(0:nb_max, 1:na_max, GIC_all(:,:,m)); set(gca,'YDir','normal');
    colorbar; xlabel('n_b'); ylabel('n_a'); 
    title(sprintf('GIC (%s)', methods{m})); hold on;
    [~, im] = min(GIC_all(:,:,m), [], 'all', 'linear');
    [na_m, nb_m] = ind2sub([na_max, nb_max+1], im);
    plot(nb_m-1, na_m, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
end
exportgraphics(fig2b, fullfile(picDir, '3_GIC_All_Methods.png'), 'Resolution', 300);

% Summary comparison plot: All Whiteness in one figure (4x2 grid)
fig3 = figure('Units', 'centimeters', 'Position', [5 5 20 18]);
for m = 1:n_methods
    subplot(4,2,m); 
    imagesc(0:nb_max, 1:na_max, whiteness_all(:,:,m)*100); set(gca,'YDir','normal');
    colorbar; xlabel('n_b'); ylabel('n_a'); 
    title(sprintf('Whiteness %% (%s)', methods{m})); hold on;
    % Mark AICc optimal (green) and GIC optimal (magenta)
    [~, im_aic] = min(AIC_all(:,:,m), [], 'all', 'linear');
    [na_aic, nb_aic] = ind2sub([na_max, nb_max+1], im_aic);
    [~, im_gic] = min(GIC_all(:,:,m), [], 'all', 'linear');
    [na_gic, nb_gic] = ind2sub([na_max, nb_max+1], im_gic);
    plot(nb_aic-1, na_aic, 'go', 'MarkerSize', 8, 'LineWidth', 2);
    plot(nb_gic-1, na_gic, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
end
exportgraphics(fig3, fullfile(picDir, '3_Whiteness_All_Methods.png'), 'Resolution', 300);

% Print summary table
fprintf('\n=== Model Order Selection Summary ===\n');
fprintf('%-6s | AICc: na,nb (AICc, W%%) | GIC: na,nb (GIC, W%%)\n', 'Method');
fprintf('-------+------------------------+----------------------\n');
for m = 1:n_methods
    [~, im_aic] = min(AIC_all(:,:,m), [], 'all', 'linear');
    [na_aic, nb_aic] = ind2sub([na_max, nb_max+1], im_aic);
    [~, im_gic] = min(GIC_all(:,:,m), [], 'all', 'linear');
    [na_gic, nb_gic] = ind2sub([na_max, nb_max+1], im_gic);
    fprintf('%-6s | %d,%d (%.1f, %.0f%%)     | %d,%d (%.1f, %.0f%%)\n', methods{m}, ...
        na_aic, nb_aic-1, AIC_all(na_aic, nb_aic, m), whiteness_all(na_aic, nb_aic, m)*100, ...
        na_gic, nb_gic-1, GIC_all(na_gic, nb_gic, m), whiteness_all(na_gic, nb_gic, m)*100);
end

% Generate LaTeX table code (AICc-based)
fprintf('\n=== LaTeX Table Code - AICc (copy to content.tex) ===\n');
fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf('\\begin{tabular}{|l|c|c|c|c|}\n');
fprintf('\\hline\n');
fprintf('\\textbf{Method} & $n_a^*$ & $n_b^*$ & \\textbf{AICc} & \\textbf{Whiteness (\\%%)} \\\\\n');
fprintf('\\hline\n');
for m = 1:n_methods
    [~, im] = min(AIC_all(:,:,m), [], 'all', 'linear');
    [na_m, nb_m] = ind2sub([na_max, nb_max+1], im);
    fprintf('%s & %d & %d & %.2f & %.1f \\\\\n', methods{m}, na_m, nb_m-1, ...
        AIC_all(na_m, nb_m, m), whiteness_all(na_m, nb_m, m)*100);
    fprintf('\\hline\n');
end
fprintf('\\end{tabular}\n');
fprintf('\\caption{Optimal model orders $(n_a^*, n_b^*)$ selected by minimizing AICc for each estimation method.}\n');
fprintf('\\label{tab:model_order_aicc}\n');
fprintf('\\end{table}\n');

% Generate LaTeX table code (GIC-based)
fprintf('\n=== LaTeX Table Code - GIC (copy to content.tex) ===\n');
fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf('\\begin{tabular}{|l|c|c|c|c|}\n');
fprintf('\\hline\n');
fprintf('\\textbf{Method} & $n_a^*$ & $n_b^*$ & \\textbf{GIC} & \\textbf{Whiteness (\\%%)} \\\\\n');
fprintf('\\hline\n');
for m = 1:n_methods
    [~, im] = min(GIC_all(:,:,m), [], 'all', 'linear');
    [na_m, nb_m] = ind2sub([na_max, nb_max+1], im);
    fprintf('%s & %d & %d & %.2f & %.1f \\\\\n', methods{m}, na_m, nb_m-1, ...
        GIC_all(na_m, nb_m, m), whiteness_all(na_m, nb_m, m)*100);
    fprintf('\\hline\n');
end
fprintf('\\end{tabular}\n');
fprintf('\\caption{Optimal model orders $(n_a^*, n_b^*)$ selected by minimizing GIC for each estimation method.}\n');
fprintf('\\label{tab:model_order_gic}\n');
fprintf('\\end{table}\n');

fprintf('\nModel order selection figures saved to %s\n', picDir);
