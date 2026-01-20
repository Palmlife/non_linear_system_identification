clear; close all; clc;
fs = 4000;
fmax = 1200;  % Max excited frequency

%% Non-parametric estimation at last power level (21)
powerLevel = 21;
[G, f_exc, varG, Y_exc, Y_even, f_even, Y_odd, f_odd, Y_noise, f_noise] = ...
    fastMethod("results_5_fast_power/outOddM", fs, powerLevel);

%% Plot 1: Output spectrum Y with bin classification
figure;
plot(f_exc, db(abs(Y_exc)), 'ko', 'MarkerSize', 4, 'LineWidth', 1.5);
hold on;
plot(f_even, db(abs(Y_even)), 'go', 'MarkerSize', 4, 'LineWidth', 1.5);
plot(f_odd, db(abs(Y_odd)), 'ro', 'MarkerSize', 4, 'LineWidth', 1.5);
plot(f_noise, db(abs(Y_noise)), 'b.', 'MarkerSize', 3);

xlabel('Frequency [Hz]');
ylabel('|Y| [dBV]');
xlim([0 fmax]);
legend('Excited (Linear)', 'Even distortions', 'Odd distortions', 'Noise', 'Location', 'best');
title(sprintf('Output Spectrum Y'));
grid on;

%% Plot 2: FRF G - Bode plot (magnitude and phase)
figure;

subplot(211);
plot(f_exc, db(abs(G)), 'b.', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Frequency [Hz]');
ylabel('|G| [dB]');
xlim([0 fmax]);
title(sprintf('FRF Estimate G'));
grid on;

subplot(212);
plot(f_exc, unwrap(angle(G))*180/pi, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Frequency [Hz]');
ylabel('Phase [°]');
xlim([0 fmax]);
grid on;
