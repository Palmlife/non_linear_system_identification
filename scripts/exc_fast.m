clear; clc; close all
% Configuration
filename = 'oddM';
outputDir = '../signals';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% --- SETTINGS TO FIX: FREQUENCY DEFINITION ---
N = 5000;               % Samples per period (Npp)
Fs = 4000;             % Sampling Frequency in Hz (to be chosen)
F_start_Hz = 1;        % Start Frequency in Hz
F_stop_Hz = 1200;       % Stop Frequency in Hz (Must be < Fs/2)
Desired_RMS = 0.1;     % Desired RMS value
R_realizations = 2;     % Number of Realizations

% 1. Convert Hz to Frequency Bins
% Fundamental frequency f0 = Fs / N
f0 = Fs / N;
f_min = round(F_start_Hz / f0);
f_max = round(F_stop_Hz / f0);

% Ensure f_min is at least 1 (DC offset bin 0 is usually avoided)
if f_min < 1; f_min = 1; end

% 2. Generate "Fast Method" Grid (Sparse Odd Multisine)
% Start with all odd bins in the calculated range 
% We verify the start bin is odd; if even, add 1.
if mod(f_min, 2) == 0; f_min = f_min + 1; end
all_odds = (f_min:2:f_max)'; 

% Randomly remove 1 out of every 4 odd lines to create detection holes [cite: 721]
group_size = 4;
Sel = [];
for i = 1:group_size:length(all_odds)
    if i + group_size - 1 <= length(all_odds)
        group = all_odds(i : i+group_size-1);
        remove_idx = randi(group_size); % Pick one to remove
        group(remove_idx) = []; 
        Sel = [Sel; group];
    else
        Sel = [Sel; all_odds(i:end)];
    end
end

% 3. Prepare Excitation Matrix (Npp x R)
Sig = zeros(N,  R_realizations); 

for r = 1:R_realizations
    % Random Phases
    phases = 2 * pi * rand(size(Sel));
    
    % Construct Frequency Domain
    X = zeros(N, 1);
    % Map selected bins to FFT indices (Bin k corresponds to Matlab index k+1)
    X(Sel + 1) = exp(1j * phases); 
    X(N - Sel + 1) = conj(X(Sel + 1)); 
    
    % Time Domain & RMS Scaling
    x_time = real(ifft(X));
    x_time = (x_time / std(x_time)) * Desired_RMS; 

    % Assign to ODD column
    Sig(:, r) = x_time;
end

% 4. Save Files
save(fullfile(outputDir, [filename '_Sig_E0_S0.mat']), 'Sig');
save(fullfile(outputDir, [filename '_Sel_E0_S0.mat']), 'Sel');
disp(['Signals saved. Freq range: ' num2str(min(Sel)*f0) ' Hz to ' num2str(max(Sel)*f0) ' Hz.']);

% --- PLOTTING ---
figure('Name', 'Fast Method Signal (Hz Scale)');
subplot(2,1,1);
plot((0:N-1)/Fs, Sig(:,1)); 
xlabel('Time (s)'); ylabel('Amplitude'); title('Time Domain'); grid on;

subplot(2,1,2);
F_spectrum = abs(fft(Sig(:,1)));
freq_axis = (0:N-1) * f0; % Frequency axis in Hz
stem(freq_axis, F_spectrum, 'Marker', 'none'); hold on;
% Mark holes
omitted = setdiff(all_odds, Sel);
if ~isempty(omitted)
    plot(omitted * f0, zeros(size(omitted)), 'rx');
end
xlim([0 F_stop_Hz * 1.1]);
xlabel('Frequency (Hz)'); ylabel('Magnitude'); 
title('Frequency Domain (Red X = Detection Lines)'); grid on;