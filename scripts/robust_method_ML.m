% Robust identification method - Version 2
% Instead of averaging FRFs, we average Y and U separately:
% 1. Average Y and U over periods
% 2. Align realizations by dividing by reference phase of U
% 3. Average Y and U over realizations
% 4. Compute FRF = Y_avg / U_avg

function [G_ML, f, noise_var, total_var, distortion_var] = robust_method_ML(file_name, fs, showPlot)
    % Load data
    [u, y, ~, sig, realizations, ~] = acquisition(file_name);

    N = size(u, 1);
    Npp = size(sig, 1); % number of samples of the original period
    periods = N / Npp; % number of periods 

    transientPeriods = 1;
    assert(transientPeriods < periods, 'Transient periods to remove exceed total number of periods.');

    u = u(transientPeriods*Npp + 1:end, :, :);
    y = y(transientPeriods*Npp + 1:end, :, :);

    % Update number of samples
    N = size(u, 1);
    periods = periods - transientPeriods;

    %% --------------- FRF estimate and variance computation ---------------

    % Preallocate
    U_avg_per_real = zeros(realizations, Npp);  % U averaged over periods
    Y_avg_per_real = zeros(realizations, Npp);  % Y averaged over periods
    G_per_real = zeros(realizations, Npp);      % G per realization
    var_Y_noise_per_real = zeros(realizations, Npp);  % Noise variance on Y per realization

    % For each realization
    for r = 1:realizations
        U_sum = zeros(Npp, 1);
        Y_sum = zeros(Npp, 1);
        
        % First pass: compute sums over periods
        Y_all_periods = zeros(periods, Npp);
        for p = 1:periods
            u_per = u((p-1)*Npp + 1:p*Npp, r, 1);
            y_per = y((p-1)*Npp + 1:p*Npp, r, 1);
            
            U = fft(u_per);
            Y = fft(y_per);
            
            U_sum = U_sum + U;
            Y_sum = Y_sum + Y;
            Y_all_periods(p, :) = Y.';
        end
        
        % Average over periods
        U_avg_per_real(r, :) = U_sum.' / periods;
        Y_avg_per_real(r, :) = Y_sum.' / periods;
        
        % G for this realization
        G_per_real(r, :) = Y_avg_per_real(r, :) ./ U_avg_per_real(r, :);
        
        % Noise variance on Y: variance across periods
        Y_mean_r = Y_avg_per_real(r, :);
        for p = 1:periods
            var_Y_noise_per_real(r, :) = var_Y_noise_per_real(r, :) + abs(Y_all_periods(p, :) - Y_mean_r).^2;
        end
        var_Y_noise_per_real(r, :) = var_Y_noise_per_real(r, :) / (periods - 1);
    end
    
    %% Align realizations by reference phase of U
    U_aligned = zeros(realizations, Npp);
    Y_aligned = zeros(realizations, Npp);
    
    for r = 1:realizations
        ref_phase = exp(1j * angle(U_avg_per_real(r, :)));
        U_aligned(r, :) = U_avg_per_real(r, :) ./ ref_phase;
        Y_aligned(r, :) = Y_avg_per_real(r, :) ./ ref_phase;
    end
    
    %% Average over realizations and compute FRF
    U_final = mean(U_aligned, 1).';
    Y_final = mean(Y_aligned, 1).';
    G_ML = Y_final ./ U_final;
    
    %% Variance estimation
    
    % 1. Noise variance on Y (averaged over realizations)
    var_Y_noise = mean(var_Y_noise_per_real, 1).';
    
    % 2. Variance on G due to noise: var_Y / |U|^2
    var_G_noise_per_real = var_Y_noise_per_real ./ (abs(U_avg_per_real).^2);
    noise_var = mean(var_G_noise_per_real, 1).' / periods;
    
    % 3. Total variance of G_BLA: variance across realizations
    var_G_total = zeros(Npp, 1);
    for r = 1:realizations
        var_G_total = var_G_total + abs(G_per_real(r, :).' - G_ML).^2;
    end

    % we're computing the variance of the mean, so divide by realizations^2
    var_G_total = var_G_total / (realizations - 1)*realizations;
    total_var = var_G_total / realizations^2;

    % 4. Stochastic nonlinear distortion variance (on Y)
    var_distortion = (total_var - noise_var) .* abs(U_final).^2;

    var_distortion(var_distortion < 0) = 0;  % Ensure non-negative
    
    % Distortion variance on G
    distortion_var = max(total_var - noise_var, 0);
    
    f = (0:Npp-1)' * (fs/Npp);

    if showPlot
        figure;
      
        subplot(211);
        set(groot, 'defaultAxesColorOrder', colororder("sail"));
        plot(f, db(G_ML), '.', 'LineWidth', 2);
        title('Robust Method v2 (Averaged Spectra)');
        xlabel('Frequency [Hz]');
        ylabel('|G| [dB]');
        xlim([1/fs fs/4]);
        ylim([min(db(G_ML))-1 1.5*max(db(G_ML))])
        grid on;
        legend('G_{BLA}')

        subplot(212);
        hold on;
        set(groot, 'defaultAxesColorOrder', colororder("sail"));
        plot(f, db(noise_var), '.','LineWidth', 2);
        plot(f, db(distortion_var),'.', 'LineWidth', 2);
        plot(f, db(total_var), '.','LineWidth', 2);
        %plot(f, db(var_Y_noise), 'LineWidth', 2);
        title('Variance Estimates');
        xlabel('Frequency [Hz]');
        ylabel('Variance [dBV^2]');
        xlim([1/fs fs/4]);
        grid on;
        legend('\sigma^2_G (noise)', '\sigma^2_G (distortion)', '\sigma^2_{G_{BLA}} (total)');
        %legend('\sigma^2_G (noise)', '\sigma^2_G (distortion)', '\sigma^2_{G_{BLA}} (total)', '\sigma^2_Y (noise)');
    end
end
