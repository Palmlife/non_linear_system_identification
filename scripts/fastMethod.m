function [G, f_exc, varG, Y_exc, Y_even, f_even, Y_odd, f_odd, Y_noise, f_noise] = fastMethod(dataFile, fs, powerLevel)
% [G, f_exc, varG, Y_exc, Y_even, f_even, Y_odd, f_odd, Y_noise, f_noise] = fastMethod(dataFile, fs, powerLevel)
%
% Computes the FRF estimate using the fast method.
% Based on a single realization, the fft bins of the output signal spectrum 
% are split in linear contribution, even nonlinear distortion, odd nonlinear 
% distortion and noise bins. 
%
% Inputs:
%   dataFile - folder name inside results/
%   fs - sampling frequency (Hz)
%   powerLevel - power level index (1-based, default: 1)
%
% Outputs:
%   G - FRF estimate at excited frequencies
%   f_exc - frequencies at which FRF is estimated
%   varG - variance estimate of the FRF (from noise floor)
%   Y_exc - output spectrum at excited frequencies
%   Y_even - output spectrum at even distortion bins
%   f_even - frequencies of even distortion bins
%   Y_odd - output spectrum at odd distortion bins  
%   f_odd - frequencies of odd distortion bins
%   Y_noise - output spectrum at noise bins
%   f_noise - frequencies of noise bins

    if nargin < 3
        powerLevel = 1;
    end

    %% Load data
    [u, y, sel, sig, ~, power_levels] = acquisition(dataFile);
    
    assert(powerLevel >= 1 && powerLevel <= power_levels, ...
        sprintf('Power level must be between 1 and %d', power_levels));

    N = size(u, 1);
    periodN = size(sig, 1); % number of samples of the original period
    repNumber = N / periodN; % number of periods in the acquired signal

    %% Remove transient
    % transient removal
        transientPeriods = 0;
        assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

        u = u(transientPeriods*periodN + 1:end, :, :);
        y = y(transientPeriods*periodN + 1:end, :, :);

    % update sizes
        N = size(y, 1);
        repNumber = N / periodN; 

    %% --------------- FFT and bin classification ---------------

    % Use second realization and selected power level
    U = fft(u(:, 2, powerLevel));
    Y = fft(y(:, 2, powerLevel));

    % bin categories
        % all excited bins, correct for number of periods
        N_exc = (sel(:, 1)*repNumber + 1)';
        % odd non linear dist: the odd bins that are not excited
        N_odd = setdiff((1:2:N/repNumber)*repNumber + 1, N_exc);
        % even non linear dist: the even bins that are not excited
        N_even = setdiff((0:2:(N-1)/repNumber)*repNumber + 1, N_exc);
        % noise bins: all other bins (due to increase in period 9)
        N_noise = setdiff(1:N, [N_exc, N_odd, N_even]);

    %% --------------- FRF Estimation ---------------
    
    % Frequency vector
    f = (0:N-1)'*(fs/N);
    
    % FRF at excited frequencies
    G = Y(N_exc) ./ U(N_exc);
    f_exc = f(N_exc);
    Y_exc = Y(N_exc);  % Output at excited frequencies
    
    % Even nonlinear distortions
    Y_even = Y(N_even);
    f_even = f(N_even);
    
    % Odd nonlinear distortions
    Y_odd = Y(N_odd);
    f_odd = f(N_odd);
    
    % Noise
    Y_noise = Y(N_noise);
    f_noise = f(N_noise);
    
    % Variance estimate from noise floor
    noise_var = var(abs(Y_noise));
    varG = noise_var ./ (abs(U(N_exc)).^2);

end