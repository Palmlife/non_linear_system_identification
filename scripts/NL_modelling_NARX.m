clear; close all; clc;

%% 1. Data Acquisition
dataFolder = "results_3_power/";
target_power_idx = 5;

[u_all, y_all, ~, sig, realizations, power_levels] = acquisition(dataFolder);

% Basic Dimensions
Npp = size(sig, 1);
N_total = size(u_all, 1);
N_periods = N_total / Npp;
N_realizations = realizations; 

fprintf("Detected %d realizations, %d periods, samples per period %d\n", ...
    N_realizations, N_periods, Npp);

%% 2. Pre-processing (Robust Method Logic)
transientPeriods = 2; 

u_raw_long = u_all(:, :, target_power_idx); 
y_raw_long = y_all(:, :, target_power_idx);

u_reshaped = reshape(u_raw_long, Npp, N_periods, N_realizations);
y_reshaped = reshape(y_raw_long, Npp, N_periods, N_realizations);

if N_periods > transientPeriods
    u_steady = u_reshaped(:, (transientPeriods+1):end, :);
    y_steady = y_reshaped(:, (transientPeriods+1):end, :);
else
    error('Not enough periods to discard transients.');
end

u_avg = squeeze(mean(u_steady, 2)); 
y_avg = squeeze(mean(y_steady, 2));

%% 3. Alignment / Synchronization (CRITICAL)
% Calculate delay lag between input and output
delay_lag = finddelay(u_avg(:,1), y_avg(:,1));
fprintf('Detected Delay: %d samples\n', delay_lag);

if delay_lag > 0
    y_avg = circshift(y_avg, -delay_lag, 1);
elseif delay_lag < 0
    u_avg = circshift(u_avg, delay_lag, 1);
end

d_check = finddelay(u_avg(:,1), y_avg(:,1));
fprintf('Delay after compensation: %d samples\n', d_check);

%% 4. NARX Parameters & Matrix Construction (No Normalization)
na = 2;          % Feedback lags
nb = 2;          % Input lags
% N_order is implicitly handled by manually adding cubic terms below

J_total = [];
y_target_total = [];

fprintf('Building NARX regressor matrix (One-Step Ahead)...\n');

for m = 1:N_realizations
    u_m = u_avg(:, m);
    y_m = y_avg(:, m);
    
    start_idx = max(na, nb) + 1;
    
    current_J = [];
    
    % --- Linear Terms ---
    % y(k-1), y(k-2)
    for i = 1:na
        current_J = [current_J, y_m(start_idx-i : end-i)];
    end
    % u(k), u(k-1), u(k-2)
    for j = 0:nb
        current_J = [current_J, u_m(start_idx-j : end-j)];
    end
    
    % --- Nonlinear Terms (Silverbox Physics) ---
    % y(k-1)^2 and y(k-1)^3
    %current_J = [current_J, y_m(start_idx-1 : end-1).^2];
    %current_J = [current_J, y_m(start_idx-1 : end-1).^3];
    
    J_total = [J_total; current_J];
    y_target_total = [y_target_total; y_m(start_idx:end)];
end

%% 5. Identification
fprintf('Solving for NARX coefficients...\n');
theta = J_total \ y_target_total;

%% 6. Validation: One-Step Ahead Prediction
% We validate on Realization #1 using the SAME regressor structure 
% (using MEASURED past outputs, not simulated ones).

u_val = u_avg(:, 2);
y_meas_val = y_avg(:, 2);
N_val = length(u_val);
start_idx = max(na, nb) + 1;

% Re-build Regressor for Validation Data
J_val = [];
% Linear Terms
for i = 1:na
    J_val = [J_val, y_meas_val(start_idx-i : end-i)];
end
for j = 0:nb
    J_val = [J_val, u_val(start_idx-j : end-j)];
end
% Nonlinear Terms
%J_val = [J_val, y_meas_val(start_idx-1 : end-1).^2];
%J_val = [J_val, y_meas_val(start_idx-1 : end-1).^3];

% Predict (OSA)
y_mod_val = J_val * theta;

% Align measurement for plotting
y_meas_plot = y_meas_val(start_idx:end);
u_val_plot  = u_val(start_idx:end);

%% 7. Visualization & Diagnostics (GMP Style)

% --- A. Calculate Metrics ---
error_sig = y_meas_plot - y_mod_val;
nmse = sum(abs(error_sig).^2) / sum(abs(y_meas_plot).^2);
nmse_db = 10*log10(nmse);
fprintf('one-step Ahead NMSe : %.4f (%.2f dB)\n', nmse, nmse_db);

% --- B. Time Domain Comparison ---
figure('Name', 'NARX OSA Results: Time Domain', 'Color', 'w');
subplot(2,1,1);
N_plot = min(1000, length(y_meas_plot));
t_idx = 1:N_plot;

plot(t_idx, y_meas_plot(t_idx), 'b', 'LineWidth', 1.5); hold on;
plot(t_idx, y_mod_val(t_idx), 'r--', 'LineWidth', 1.2);
plot(t_idx, error_sig(t_idx), 'k', 'LineWidth', 0.5);
legend('Measurement', 'NARX (OSA)', 'Error');
title(sprintf('Time Domain Fit (First %d samples)', N_plot));
grid on; ylabel('Amplitude');

subplot(2,1,2);
plot(y_meas_plot, 'b'); hold on;
plot(y_mod_val, 'r--');
title(sprintf('Full Sequence Comparison (NMSE: %.2f dB)', nmse_db));
grid on; xlabel('Sample'); ylabel('Amplitude');

% --- C. Frequency Domain (PSD) ---
L = length(y_meas_plot);
window = hann(L); 
Y_meas_dB = 10*log10(fftshift(abs(fft(y_meas_plot .* window)/L).^2));
Y_mod_dB  = 10*log10(fftshift(abs(fft(y_mod_val  .* window)/L).^2));
E_dB      = 10*log10(fftshift(abs(fft(error_sig .* window)/L).^2));
f_axis    = linspace(-0.5, 0.5, L);

figure('Name', 'NARX OSA Results: Frequency Domain', 'Color', 'w');
plot(f_axis, Y_meas_dB, 'b', 'LineWidth', 1.5); hold on;
plot(f_axis, Y_mod_dB, 'r--', 'LineWidth', 1.2);
plot(f_axis, E_dB, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
legend('Measurement', 'Model', 'Error Floor');
title('Power Spectral Density');
xlabel('Normalized Frequency'); ylabel('Power (dB)');
grid on; ylim([-120 0]);

% --- D. Hysteresis / AM-AM Plot ---
figure('Name', 'NARX OSA Results: Hysteresis (AM-AM)', 'Color', 'w');

step = max(1, floor(length(u_val_plot)/5000)); 
idx_scatter = 1:step:length(u_val_plot);

scatter(u_val_plot(idx_scatter), y_meas_plot(idx_scatter), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3); hold on;
scatter(u_val_plot(idx_scatter), y_mod_val(idx_scatter),  10, 'r', 'filled', 'MarkerFaceAlpha', 0.3);

legend('Measurement', 'NARX (OSA)', 'Location', 'Best');
xlabel('Input u(t) (Force)');
ylabel('Output y(t) (Displacement)');
title('System Hysteresis Loop (One-Step Ahead)');
grid on;