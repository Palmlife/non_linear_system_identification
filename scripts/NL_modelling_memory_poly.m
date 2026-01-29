clear; close all; clc;

%% 1. Data Acquisition
dataFolder = "results_3_power/";
target_power_idx = 5; % Example: Choose 5th power level

% Load data
[u_all, y_all, ~, sig, realizations, power_levels] = acquisition(dataFolder);

% Basic Dimensions
Npp = size(sig, 1);           % Samples per period
N_total = size(u_all, 1);     % Total samples (all periods)
N_periods = N_total / Npp;    % Number of periods
N_realizations = realizations; 

fprintf("Detected %d realizations, %d periods, samples per period %d\n", ...
    N_realizations, N_periods, Npp);

%% 2. User Parameters
N_order = 3;      % Nonlinearity order (e.g., 1, 3, 5...)
M_lag = 100;        % Memory depth
transientPeriods = 2; % Number of initial periods to discard


%% 3. Pre-processing (Robust Method Logic)
% Extract data for the specific power level
% u_all is (Samples, Realizations, PowerLevels)
u_raw_long = u_all(:, :, target_power_idx); 
y_raw_long = y_all(:, :, target_power_idx);

% Reshape to (Npp, Periods, Realizations)
u_reshaped = reshape(u_raw_long, Npp, N_periods, N_realizations);
y_reshaped = reshape(y_raw_long, Npp, N_periods, N_realizations);

% Discard Transient Periods (e.g., period 1)
% Keep indices from (transientPeriods + 1) to end
if N_periods > transientPeriods
    u_steady = u_reshaped(:, (transientPeriods+1):end, :);
    y_steady = y_reshaped(:, (transientPeriods+1):end, :);
else
    error('Not enough periods to discard transients.');
end

% Average over the remaining periods (Robust Method) to remove noise
% Result is (Npp, 1, Realizations) -> Squeeze to (Npp, Realizations)
u_avg = squeeze(mean(u_steady, 2)); 
y_avg = squeeze(mean(y_steady, 2));

% u_avg is now clean, steady-state data with one period per realization.

%% INSERT THIS BLOCK AFTER AVERAGING AND BEFORE BUILDING J

% 1. Calculate the delay lag between input and output
delay_lag = finddelay(u_avg(:,1), y_avg(:,1));

fprintf('Detected Delay: %d samples\n', delay_lag);

% 2. Align the signals
if delay_lag > 0
    % Output is delayed -> Shift it BACK (remove the start, pad the end)
    % Or simply circular shift if PISPO (Periodic)
    y_avg = circshift(y_avg, -delay_lag, 1);
elseif delay_lag < 0
    % Input is delayed (rare) -> Shift input
    u_avg = circshift(u_avg, delay_lag, 1);
end

% 3. Verify Alignment (Optional but recommended)
% The cross-correlation peak should now be at lag 0
d_check = finddelay(u_avg(:,1), y_avg(:,1));
fprintf('Delay after compensation: %d samples (Should be 0)\n', d_check);

%% CONTINUE TO STEP 4 (Building J) ...


%% 4. Construct Regressor Matrix (J)
J_total = [];
y_total = [];

fprintf('Building regressor matrix for %d realizations...\n', N_realizations);

for m = 1:N_realizations
    % Extract the averaged period for this realization
    u_m = u_avg(:, m);
    y_m = y_avg(:, m);
    
    % Build MP matrix for this single period
    J_m = build_MP_matrix(u_m, N_order, M_lag);
    
    % IMPORTANT: Remove the first M_lag samples from the *regression* % because the history u(m-l) wraps around or is zero-padded incorrectly 
    % at the start of the vector.
    valid_indices = (M_lag + 1):Npp;
    
    J_total = [J_total; J_m(valid_indices, :)];
    y_total = [y_total; y_m(valid_indices)];
end

%% 5. Identification (Linear Least Squares)
% Solve y = J * A
fprintf('Solving for coefficients...\n');
A = J_total \ y_total; 

%% 6. Validation on First Realization
u_val = u_avg(:, 1);
y_val_meas = y_avg(:, 1);

J_val = build_MP_matrix(u_val, N_order, M_lag);
y_val_model = J_val * A;

% Apply the same mask to validation to ensure fair comparison
mask = (M_lag + 1):Npp;
y_val_meas_clipped = y_val_meas(mask);
y_val_model_clipped = y_val_model(mask);

nmse_db = 10*log10(sum(abs(y_val_meas_clipped - y_val_model_clipped).^2) / ...
                   sum(abs(y_val_meas_clipped).^2));

fprintf('NMSE: %.2f dB\n', nmse_db);

figure;
plot(abs(y_val_meas_clipped)); hold on;
plot(abs(y_val_model_clipped), '--');
legend('Measurement', 'Model');
title(['Model Fit (NMSE: ' num2str(nmse_db) ' dB)']);

%% --- Diagnostic Plots for Goodness of Fit ---

% 1. Frequency Domain Comparison (PSD)
% ------------------------------------
L = length(y_val_meas);
window = hann(L); % Use windowing for cleaner spectrum
Y_meas_dB = 10*log10(fftshift(abs(fft(y_val_meas .* window)/L).^2));
Y_mod_dB  = 10*log10(fftshift(abs(fft(y_val_model .* window)/L).^2));
Error_dB  = 10*log10(fftshift(abs(fft((y_val_meas - y_val_model) .* window)/L).^2));
freq_axis = linspace(-0.5, 0.5, L); % Normalized frequency

figure('Name', 'Frequency Domain Analysis');
plot(freq_axis, Y_meas_dB, 'LineWidth', 1.2); hold on;
plot(freq_axis, Y_mod_dB, '--', 'LineWidth', 1.2);
plot(freq_axis, Error_dB, 'Color', [0.5 0.5 0.5]); % Grey for error
legend('Measurement', 'Model', 'Error PSD');
grid on;
title('Output Power Spectral Density');
xlabel('Normalized Frequency'); ylabel('Power (dB)');
ylim([-100 0]); 

% 2. AM-AM and AM-PM Analysis
% ---------------------------
% Calculate magnitudes and phases
u_mag = abs(u_val);
y_meas_mag = abs(y_val_meas);
y_mod_mag  = abs(y_val_model);

% Phase difference (unwrapped) in degrees
% Note: Simple subtraction works if freq is baseband. 
% If signals are complex envelopes, this is valid.
phase_meas = rad2deg(angle(y_val_meas) - angle(u_val));
phase_mod  = rad2deg(angle(y_val_model) - angle(u_val));

% Wrap phase to [-180, 180] for cleaner plotting if needed
phase_meas = mod(phase_meas + 180, 360) - 180;
phase_mod  = mod(phase_mod  + 180, 360) - 180;

figure('Name', 'AM-AM and AM-PM');

% AM-AM Plot
subplot(1,2,1);
plot(u_mag, y_meas_mag, '.', 'MarkerSize', 1, 'Color', [0 0.4470 0.7410]); hold on;
plot(u_mag, y_mod_mag, '.', 'MarkerSize', 1, 'Color', [0.8500 0.3250 0.0980]);
xlabel('|u(t)| (Input Amplitude)');
ylabel('|y(t)| (Output Amplitude)');
title('AM-AM Characteristic');
legend('Measured', 'Modelled', 'Location', 'best');
grid on;

% AM-PM Plot
subplot(1,2,2);
plot(u_mag, phase_meas, '.', 'MarkerSize', 1, 'Color', [0 0.4470 0.7410]); hold on;
plot(u_mag, phase_mod, '.', 'MarkerSize', 1, 'Color', [0.8500 0.3250 0.0980]);
xlabel('|u(t)| (Input Amplitude)');
ylabel('Phase Change (deg)');
title('AM-PM Characteristic');
ylim([-20 20]); % Adjust based on your actual phase variation
grid on;


%% ---------------------------------------------------------
%  Local Function: Memory Polynomial Matrix Builder
%  ---------------------------------------------------------
function J = build_MP_matrix(u, N_order, M_lag)
    N = length(u);
    % Total coefficients = (M_lag + 1) * N_order
    J = zeros(N, (M_lag + 1) * N_order);
    
    col_idx = 1;
    
    for l = 0:M_lag
        % Shift input by lag l (fill with 0 at start)
        u_lagged = [zeros(l, 1); u(1:end-l)];
        u_abs = abs(u_lagged);
        
        for n = 1:N_order
            % Basis function: u(m-l) * |u(m-l)|^(n-1)
            if n == 1
                term = u_lagged;
            else
                term = u_lagged .* (u_abs .^ (n - 1));
            end
            
            J(:, col_idx) = term;
            col_idx = col_idx + 1;
        end
    end
end