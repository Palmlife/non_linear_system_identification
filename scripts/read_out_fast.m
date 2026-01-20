clear; clc;
%% ------- Output display (Fast Method) ----------
Fs = 4000;

[u, y, sel, sig, ~, ~] = acquisition("results_5_fast_power/", "outOddM");

N = size(u, 1);
periodN = size(sig, 1); % number of samples of the original period
repNumber = N / periodN; % number of periods in the acquired signal

%% Remove transient
transientPeriods = 1;
assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

u = u(transientPeriods*periodN + 1:end, :, :);
y = y(transientPeriods*periodN + 1:end, :, :);

% update sizes
N = size(y, 1);
repNumber = N / periodN;

%% FFT
U = fft(u(:, 1));
Y = fft(y(:, 1));

%% Bin categories (fast method: odd multisine excitation)
N_exc = (sel(:, 1)*repNumber + 1)';
N_odd = setdiff((1:2:N/repNumber)*repNumber + 1, N_exc);
N_even = setdiff((0:2:(N-1)/repNumber)*repNumber + 1, N_exc);
N_noise = setdiff(1:N, [N_exc, N_odd, N_even]);

%% Plots
f = (0:N-1)'*(Fs/N);

figure;
sgtitle('Fast Method - Output Spectrum Analysis');

% Input spectrum
subplot(211);
plot(f, db(abs(U)), '.', 'LineWidth', 0.5);
xlim([0 Fs/2]);
title('Input spectrum');
xlabel('Frequency [Hz]');
ylabel('Amplitude [dBV]');
grid on;

% Output spectrum with bin classification
subplot(212);
indices = {N_exc, N_even, N_odd, N_noise};
colors = {'k', 'g', 'r', 'b'};
labels = {'Excited (Linear)', 'Even distortions', 'Odd distortions', 'Noise'};

for itr = 1:4
    plot(f(indices{itr}), db(abs(Y(indices{itr}))), 'o', 'Color', colors{itr}, 'MarkerSize', 4);
    hold on;
end

xlim([0 Fs/2]);
title('Output spectrum (classified bins)');
xlabel('Frequency [Hz]');
ylabel('Amplitude [dBV]');
legend(labels, 'Location', 'best');
grid on;
