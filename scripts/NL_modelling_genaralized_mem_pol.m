clear; close all; clc;

%% 1. Data Acquisition (Standard)
fprintf('Acquiring data...\n');

dataFolder = "results_3_power/";
target_power_idx = 10; 

[u_all, y_all, ~, sig, realizations, power_levels] = acquisition(dataFolder);

% Basic Dimensions
Npp = size(sig, 1);
N_total = size(u_all, 1);
N_periods = N_total / Npp;
N_realizations = realizations; 

% Pre-processing: Transient Removal & Averaging
transientPeriods = 1;
u_raw = reshape(u_all(:, :, target_power_idx), Npp, N_periods, N_realizations);
y_raw = reshape(y_all(:, :, target_power_idx), Npp, N_periods, N_realizations);

if N_periods > transientPeriods
    u_avg = squeeze(mean(u_raw(:, (transientPeriods+1):end, :), 2));
    y_avg = squeeze(mean(y_raw(:, (transientPeriods+1):end, :), 2));
else
    error('Not enough periods.');
end

%% 2. User Parameters for GMP (Equation 3.72)
% Aligned (Standard MP) part
N_a = 3;      % Nonlinear order for aligned terms
M_a = 100;     % Memory depth for aligned terms

% Lagging Cross-terms (Signal leads Envelope)
N_b = 4;      % Nonlinear order
M_b = 100;      % Memory depth
V_lag = 10;    % Max lag shift (v)

% Leading Cross-terms (Signal lags Envelope)
N_c = 4;      % Nonlinear order
M_c = 100;      % Memory depth
W_lead = 10;   % Max lead shift (w)

%% 3. Construct GMP Regressor Matrix (Equation 3.72)
J_total = [];
y_total = [];

fprintf('Building GMP regressor matrix...\n');

for m = 1:N_realizations
    u_m = u_avg(:, m);
    y_m = y_avg(:, m);
    
    % We need to respect the maximum history required by any term
    max_lag = max([M_a, M_b + V_lag, M_c]); 
    % Note: Leading terms look *forward* in the envelope buffer, 
    % but physically we usually implement "leading" relative to a delayed center.
    % For causal implementation, we just ensure we have enough buffer.
    
    % Standard MP Matrix (Part 1 of Eq 3.72) 
    J_aligned = build_MP_basis(u_m, N_a, M_a);
    
    % Lagging Cross-Term Matrix (Part 2 of Eq 3.72) 
    J_lagging = build_lagging_basis(u_m, N_b, M_b, V_lag);
    
    % Leading Cross-Term Matrix (Part 3 of Eq 3.72) 
    J_leading = build_leading_basis(u_m, N_c, M_c, W_lead);
    
    % Concatenate all bases
    J_m = [J_aligned, J_lagging, J_leading];
    
    % Mask invalid indices (start of signal)
    % We need to discard enough samples so that (m-l-v) is valid
    start_idx = max_lag + 1;
    valid_indices = start_idx:Npp;
    
    J_total = [J_total; J_m(valid_indices, :)];
    y_total = [y_total; y_m(valid_indices)];
end

%% 4. Identification
fprintf('Solving for coefficients (Size: %d x %d)...\n', size(J_total));
coeffs = J_total \ y_total;

%% 5. Validation (Define y_val_meas and y_val_model)

% Select the first realization for validation
u_val = u_avg(:, 1);
y_val_meas = y_avg(:, 1);

% Re-build the GMP Matrix for this single validation realization
% (Using the same N_a, M_a, N_b, M_b, etc. defined in Step 2)

% 1. Standard Aligned Part
J_aligned_val = build_MP_basis(u_val, N_a, M_a);

% 2. Lagging Cross-Terms
J_lagging_val = build_lagging_basis(u_val, N_b, M_b, V_lag);

% 3. Leading Cross-Terms
J_leading_val = build_leading_basis(u_val, N_c, M_c, W_lead);

% Concatenate
J_val = [J_aligned_val, J_lagging_val, J_leading_val];

% Apply the same valid indices mask used during identification
% (or just ensure sizes match if you want to predict the full valid range)
max_lag = max([M_a, M_b + V_lag, M_c]);
valid_indices_val = (max_lag + 1):length(u_val);

% Truncate to valid range for fair comparison
y_val_meas = y_val_meas(valid_indices_val);
u_val      = u_val(valid_indices_val);
J_val      = J_val(valid_indices_val, :);

% Calculate Model Output
y_val_model = J_val * coeffs;

%% ---------------------------------------------------------
%  Local Functions for Basis Construction 
%  ---------------------------------------------------------

% 1. Standard Aligned Terms: u(m-l) * |u(m-l)|^(n-1)
function J = build_MP_basis(u, N_order, M_depth)
    N = length(u);
    J = [];
    for l = 0:M_depth
        u_del = [zeros(l,1); u(1:end-l)];
        for n = 1:N_order
            J = [J, u_del .* (abs(u_del).^(n-1))];
        end
    end
end

% 2. Lagging Cross-Terms: u(m-l) * |u(m-l-v)|^(n-1)
function J = build_lagging_basis(u, N_order, M_depth, V_max)
    N = length(u);
    J = [];
    % Sum over l (memory)
    for l = 0:M_depth
        u_del_l = [zeros(l,1); u(1:end-l)]; % u(m-l)
        
        % Sum over v (lag depth) - starting from 1 
        for v = 1:V_max
            total_delay = l + v;
            u_del_env = [zeros(total_delay,1); u(1:end-total_delay)]; % u(m-l-v)
            
            % Sum over n (order) - starting from 2 
            for n = 2:N_order
                term = u_del_l .* (abs(u_del_env).^(n-1));
                J = [J, term];
            end
        end
    end
end

% 3. Leading Cross-Terms: u(m-l) * |u(m-l+w)|^(n-1)
function J = build_leading_basis(u, N_order, M_depth, W_max)
    N = length(u);
    J = [];
    
    % NOTE: "Leading" implies looking into the future relative to 'l'.
    % To keep this causal/implementable, 'l' usually starts large enough,
    % or we treat 'u(m-l)' as the reference and ensure u(m-l+w) exists.
    
    for l = 0:M_depth
        u_del_l = [zeros(l,1); u(1:end-l)]; % u(m-l)
        
        % Sum over w (lead depth) - starting from 1 
        for w = 1:W_max
            % To get u(m-l+w), we shift u(m-l) BACK by w (to the left)
            % Or simply: u(m - (l-w))
            eff_lag = l - w;
            
            if eff_lag >= 0
                u_del_env = [zeros(eff_lag,1); u(1:end-eff_lag)];
            else
                % Negative lag means future relative to index 1.
                % In offline processing, we can grab u(1-eff_lag : end) padded at end
                shift = -eff_lag;
                u_del_env = [u(shift+1:end); zeros(shift,1)];
            end
            
            % Sum over n (order) - starting from 2 
            for n = 2:N_order
                term = u_del_l .* (abs(u_del_env).^(n-1));
                J = [J, term];
            end
        end
    end
end

%% 6. Visualization & Plotting

% --- A. Calculate Performance Metrics ---
% Ensure vectors are column vectors
y_meas = y_val_meas(:);
y_mod  = y_val_model(:);
u_in   = u_val(:);

%conditioning check
J_cond = cond(J_val);
fprintf('Condition number of the GMP regressor matrix: %.2e\n', J_cond);

% MSE and Normalized Mean Square Error (NMSE)
error_sig = y_meas - y_mod;
mse = mean(abs(error_sig).^2);
fprintf('Mean Square Error (MSE): %.4e\n', mse);

mse_db = 10*log10(mse);
fprintf('Mean Square Error (MSE) in dB: %.2f dB\n', mse_db);

nmse = sum(abs(error_sig).^2) / sum(abs(y_meas).^2);

nmse_db = 10*log10(nmse);
fprintf('Final NMSE: %.2f dB\n', nmse_db);

% --- B. Time Domain Comparison ---
figure('Name', 'GMP Results: Time Domain', 'Color', 'w');
subplot(2,1,1);
% Plot a segment (e.g., first 500 samples) to see detail
N_plot = min(1000, length(y_meas));
t_idx = 1:N_plot;

plot(t_idx, y_meas(t_idx), 'b', 'LineWidth', 1.5); hold on;
plot(t_idx, y_mod(t_idx), 'r--', 'LineWidth', 1.2);
plot(t_idx, error_sig(t_idx), 'k', 'LineWidth', 0.5);
legend('Measurement', 'GMP Model', 'Error');
title(sprintf('Time Domain Fit (First %d samples)', N_plot));
grid on; ylabel('Amplitude');

subplot(2,1,2);
% Plot the whole sequence to check for stability
plot(y_meas, 'b'); hold on;
plot(y_mod, 'r--');
title(sprintf('Full Sequence Comparison (NMSE: %.2f dB)', nmse_db));
grid on; xlabel('Sample'); ylabel('Amplitude');


% --- C. Frequency Domain (PSD) ---
L = length(y_meas);
fs = 4000; % Sample frequency in Hz
%window = hann(L);
% Compute PSD (using windowing for dynamic range)
Y_meas = fft(y_meas);
Y_mod  = fft(y_mod);
E      = fft(error_sig);
%Y_meas = fft(y_meas .* window, L);
%Y_mod  = fft(y_mod  .* window, L);
%E      = fft(error_sig .* window, L);

f_axis = fs * (-0.5:1/L:0.5-1/L); % Frequency axis in Hz (centered)
Y_meas_dB = 10*log10(fftshift(abs(Y_meas/L).^2));
Y_mod_dB  = 10*log10(fftshift(abs(Y_mod/L).^2));
E_dB      = 10*log10(fftshift(abs(E/L).^2));

figure('Name', 'GMP Results: Frequency Domain', 'Color', 'w');
plot(f_axis, Y_meas_dB, 'b', 'LineWidth', 1.5); hold on;
plot(f_axis, Y_mod_dB, 'r--', 'LineWidth', 1.2);
plot(f_axis, E_dB, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
legend('Measurement', 'Model', 'Error Floor');
title('Power Spectral Density');
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
grid on; ylim([-120 0]);


% --- D. Hysteresis / AM-AM Plot ---
% For Silverbox: This plots Displacement (Output) vs Force (Input)
% This visualizes the nonlinear stiffness.
figure('Name', 'GMP Results: Hysteresis (AM-AM)', 'Color', 'w');

% Downsample for plotting speed if data is huge
step = max(1, floor(length(u_in)/5000)); 
idx_scatter = 1:step:length(u_in);

scatter(u_in(idx_scatter), y_meas(idx_scatter), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.3); hold on;
scatter(u_in(idx_scatter), y_mod(idx_scatter),  10, 'r', 'filled', 'MarkerFaceAlpha', 0.3);

legend('Measurement', 'GMP Model', 'Location', 'Best');
xlabel('Input u(t) (Force)');
ylabel('Output y(t) (Displacement)');
title('System Hysteresis Loop (Static Nonlinearity + Memory)');
grid on;

% Create a "difference" plot to see where the model fails in the phase plane
% (Optional: useful for mechanical systems)
figure('Name', 'Modelling Error vs Input', 'Color', 'w');
scatter(u_in(idx_scatter), error_sig(idx_scatter), 10, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Input u(t)'); ylabel('Error e(t)');
title('Residual Error Distribution');
grid on;


