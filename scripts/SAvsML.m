% Compare robust_method_SA vs robust_method_ML
clear; close all; clc;
fs = 4000; fmax = 1200;

%% Non-parametric estimation
[G1, f, nv1, tv1, dv1] = robust_method_SA("results_2/fullMOut", fs, 0);
[G2, ~, nv2, tv2, dv2] = robust_method_ML("results_2/fullMOut", fs, 0);
idx = (f > 0) & (f < fmax);

%% Compare FRF and variances
figure;
subplot(3,1,1);
plot(f, db(G1), 'b.', f, db(G2), 'r.', 'MarkerSize', 3);
xlabel('f [Hz]'); ylabel('|G| [dB]'); xlim([0 fmax]); grid on;
legend('G_{BLA} (avg FRF)', 'G_{BLA} (avg spectra)'); title('FRF Comparison');

subplot(3,1,2); hold on;
plot(f(idx), db(tv1(idx)), 'bo', 'MarkerSize', 3);
plot(f(idx), db(nv1(idx)), 'g^', 'MarkerSize', 3);
plot(f(idx), db(dv1(idx)), 'rs', 'MarkerSize', 3);
xlabel('f [Hz]'); ylabel('Variance [dB]'); xlim([0 fmax]); grid on;
legend('\sigma^2_{G} total', '\sigma^2_{G} noise', '\sigma^2_{G} distortion');
title('Variance (avg FRF)');
ax2 = gca;

subplot(3,1,3); hold on;
plot(f(idx), db(tv2(idx)), 'bo', 'MarkerSize', 3);
plot(f(idx), db(nv2(idx)), 'g^', 'MarkerSize', 3);
plot(f(idx), db(dv2(idx)), 'rs', 'MarkerSize', 3);
xlabel('f [Hz]'); ylabel('Variance [dB]'); xlim([0 fmax]); grid on;
legend('\sigma^2_{G} total', '\sigma^2_{G} noise', '\sigma^2_{G} distortion');
title('Variance (avg spectra)');
ax3 = gca;

linkaxes([ax2, ax3], 'y');

fprintf('RMS FRF difference: %.4e\n', rms(G1(idx) - G2(idx)));
