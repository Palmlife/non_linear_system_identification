clear; close all; clc;
fs = 4000;
fmax = 1200;  % Max excited frequency

%% Select 6 power levels spread across the range (1 to 21)
powerLevels = [1, 5, 9, 13, 17, 21];

%% Figure 1: Output spectrum with bin classification (6 subplots)
figure;
sgtitle('Output Spectrum - Bin Classification');

for idx = 1:6
    powerLevel = powerLevels(idx);
    
    [G, f_exc, varG, Y_even, f_even, Y_odd, f_odd, Y_noise, f_noise] = ...
        fastMethod("results_5_fast_power/outOddM", fs, powerLevel);

    subplot(2, 3, idx);
    plot(f_exc, db(abs(G .* 1)), 'ko', 'MarkerSize', 3);
    hold on;
    plot(f_even, db(abs(Y_even)), 'go', 'MarkerSize', 3);
    plot(f_odd, db(abs(Y_odd)), 'ro', 'MarkerSize', 3);
    plot(f_noise, db(abs(Y_noise)), 'b.', 'MarkerSize', 2);

    xlabel('Frequency [Hz]');
    ylabel('|Y| [dBV]');
    xlim([0 fs/2]);
    title(sprintf('Power Level %d', powerLevel));
    grid on;
end
legend('Excited', 'Even', 'Odd', 'Noise', 'Location', 'best');

%% Figure 2: FRF magnitude (6 subplots)
figure;
sgtitle('FRF Magnitude - Fast Method');

for idx = 1:6
    powerLevel = powerLevels(idx);
    
    [G, f_exc, ~, ~, ~, ~, ~, ~, ~] = ...
        fastMethod("results_5_fast_power/outOddM", fs, powerLevel);

    subplot(2, 3, idx);
    plot(f_exc, db(abs(G)), 'b.-', 'LineWidth', 1, 'MarkerSize', 4);
    xlabel('Frequency [Hz]');
    ylabel('|G| [dB]');
    xlim([0 fmax]);
    title(sprintf('Power Level %d', powerLevel));
    grid on;
end

%% Figure 3: FRF phase (6 subplots)
figure;
sgtitle('FRF Phase - Fast Method');

for idx = 1:6
    powerLevel = powerLevels(idx);
    
    [G, f_exc, ~, ~, ~, ~, ~, ~, ~] = ...
        fastMethod("results_5_fast_power/outOddM", fs, powerLevel);

    subplot(2, 3, idx);
    plot(f_exc, unwrap(angle(G))*180/pi, 'b.-', 'LineWidth', 1, 'MarkerSize', 4);
    xlabel('Frequency [Hz]');
    ylabel('Phase [°]');
    xlim([0 fmax]);
    title(sprintf('Power Level %d', powerLevel));
    grid on;
end