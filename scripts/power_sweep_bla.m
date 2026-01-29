
clear; close all; clc;
fs = 4000;
fmax = 1200;
dataFolder = "results_3_power/";

[~, ~, ~, ~, realizations, power_levels] = acquisition(dataFolder);
fprintf("Detected %d realizations and %d power levels in %s\n", realizations, power_levels, dataFolder);

%calculate the rms of the input signal at each power level

[u_all, ~, ~, sig, ~, ~] = acquisition(dataFolder);
Npp = size(sig, 1); % samples per period
N = size(u_all, 1); % total samples
periods = N / Npp;
transientPeriods = 1; % or whatever you use in robust_method_ML

% Remove transient period
u_all = u_all(transientPeriods*Npp+1:end, :, :);

% Now recalculate input_rms
input_rms = zeros(1, power_levels-1);
for p = 2:power_levels
    u_p = u_all(:, :, p); % samples x realizations
    input_rms(p-1) = sqrt(mean(u_p(:).^2));
end

%input_rms_dbm = 20*log10(input_rms/1e-3); % dBm, assuming 1 Ohm reference




% Preallocate using first run
[G_ref, f_ref, noise_ref, total_ref, distortion_ref] = robust_method_ML(dataFolder, fs, 0, 1);
K = numel(f_ref);
G_all = zeros(K, power_levels);
noise_all = zeros(K, power_levels);
dist_all = zeros(K, power_levels);
total_all = zeros(K, power_levels);

G_all(:, 1) = G_ref;
noise_all(:, 1) = noise_ref;
dist_all(:, 1) = distortion_ref;
total_all(:, 1) = total_ref;

for p = 2:power_levels
    [G_p, f_p, noise_p, total_p, distortion_p] = robust_method_ML(dataFolder, fs, 0, p);
    G_all(:, p) = G_p;
    noise_all(:, p) = noise_p;
    dist_all(:, p) = distortion_p;
    total_all(:, p) = total_p;
   
end

idx = (f_ref > 0) & (f_ref < fmax);
% Exclude the first power index (often reference/zero input)
power_idx = 1:power_levels-1;

mean_noise = mean(noise_all(idx, 2:end), 1);
mean_dist = mean(dist_all(idx, 2:end), 1);
mean_total = mean(total_all(idx, 2:end), 1);
mean_mag = mean(db(abs(G_all(idx, 2:end))), 1);

dist_ratio = mean_dist ./ max(mean_total, eps);
noise_ratio = mean_noise ./ max(mean_total, eps);

scriptDir = fileparts(mfilename("fullpath"));
picDir = fullfile(scriptDir, "..", "report", "pic");
if ~exist(picDir, "dir"); mkdir(picDir); end

cmap = turbo(power_levels);
fig1 = figure("Units", "centimeters", "Position", [4 4 16 10]);
xlabel("Frequency [Hz]");
ylabel("|G_{BLA}| [dB]");
title("BLA magnitude vs input RMS power [V^2]");
hold on;
for p = 2:power_levels
    plot(f_ref(idx), db(abs(G_all(idx, p))), "Color", cmap(p, :));
end
colormap(cmap);
legend_entries = cell(power_levels-1, 1);
for p = 2:power_levels
    legend_entries{p-1} = sprintf("Input Power: %.3f V^2", input_rms(p-1));
end
legend(legend_entries, "Location", "best");
xlim([0 200]);
exportgraphics(fig1, fullfile(picDir, "5_BLA_Magnitude_vs_Freq_and_Power_1.png"), "Resolution", 300);
hold off;

fig1 = figure("Units", "centimeters", "Position", [4 4 16 10]);
xlabel("Frequency [Hz]");
ylabel("|G_{BLA}| [dB]");
title("BLA magnitude vs input RMS power [V^2]");
hold on;
for p = 2:power_levels
    plot(f_ref(idx), db(abs(G_all(idx, p))), "Color", cmap(p, :));
end
colormap(cmap);
legend_entries = cell(power_levels-1, 1);
for p = 2:power_levels
    legend_entries{p-1} = sprintf("Input Power: %.3f V^2", input_rms(p-1));
end
legend(legend_entries, "Location", "best");
xlim([55 85]);
exportgraphics(fig1, fullfile(picDir, "5_BLA_Magnitude_vs_Freq_and_Power_2.png"), "Resolution", 300);
hold off;

% Per-frequency nonlinear distortion vs power (line overlay)
fig2 = figure("Units", "centimeters", "Position", [4 4 16 10]);
hold on;
for p = 2:power_levels
    plot(f_ref(idx), db(dist_all(idx, p)), "Color", cmap(p, :));
end
colormap(cmap);
title("Nonlinear distortion vs frequency and input power [V^2]");
% cb = colorbar;
% cb.Label.String = "Power index";
% cb.Ticks = linspace(0, 1, min(power_levels, 6));
% cb.TickLabels = round(linspace(0, power_levels-1, numel(cb.Ticks)));
xlabel("Frequency [Hz]");
ylabel("Distortion variance [dB]");
xlim([0 200]);
legend(legend_entries, "Location", "best");
exportgraphics(fig2, fullfile(picDir, "5_Distortion_vs_Freq_and_Power_1.png"), "Resolution", 300);
hold off;

% Per-frequency nonlinear distortion and noise vs power (line overlay)
fig2 = figure("Units", "centimeters", "Position", [4 4 16 10]);
hold on;
for p = 2:power_levels
   % plot(f_ref(idx), db(dist_all(idx, p)), "Color", cmap(p, :));
    plot(f_ref(idx), db(noise_all(idx, p)), "Color", cmap(p, :));
end
colormap(cmap);
title("Nonlinear distortion and noise vs frequency and input power [V^2]");
% cb = colorbar;
% cb.Label.String = "Power index";
% cb.Ticks = linspace(0, 1, min(power_levels, 6));
% cb.TickLabels = round(linspace(0, power_levels-1, numel(cb.Ticks)));
xlabel("Frequency [Hz]");
ylabel("Distortion variance [dB]");
xlim([0 200]);
legend(legend_entries, "Location", "best");
exportgraphics(fig2, fullfile(picDir, "5_Distortion_vs_Freq_and_Power_2.png"), "Resolution", 300);
hold off;

% % Heatmap view of nonlinear distortion (dB)
% fig3 = figure("Units", "centimeters", "Position", [4 4 16 10]);
% hold on;
% imagesc(input_rms, f_ref(idx), db(dist_all(idx, 2:end)));
% set(gca, "YDir", "normal");
% xlabel("Input RMS power [V^2]");
% ylabel("Frequency [Hz]");
% title("Nonlinear distortion (dB) heatmap");
% colormap(turbo);
% grid on;
% cb2 = colorbar; cb2.Label.String = "Distortion variance [dB]";
% %ylim([0 fmax]);
% exportgraphics(fig3, fullfile(picDir, "5_Distortion_Heatmap.png"), "Resolution", 300);

% hold off;


% Plot noise power vs frequency for all power levels (except the first)
fig_noise = figure("Units", "centimeters", "Position", [4 4 16 9]);
hold on;
for p = 2:power_levels
    plot(f_ref(idx), db(noise_all(idx, p)), "Color", cmap(p, :));
end
colormap(cmap);
xlabel("Frequency [Hz]");
ylabel("Noise Power [dB]");
title("Noise Power vs Frequency for Different Input Powers");
legend_entries_noise = cell(power_levels-1, 1);
for p = 2:power_levels
    legend_entries_noise{p-1} = sprintf("Input Power: %.3f V^2", input_rms(p-1));
end
%legend(legend_entries_noise);
xlim([0 200]);
grid on;
exportgraphics(fig_noise, fullfile(picDir, "5_Noise_vs_Freq_and_Power.png"), "Resolution", 300);
hold off;


fprintf("Figures saved to %s\\n", picDir);
