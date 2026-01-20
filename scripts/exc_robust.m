% Configuration
filename = 'fullM';
outputDir = '../signals';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% --- SETTINGS: FREQUENCY & RMS ---
N = 5000;               % Samples per period (Npp)
Fs = 4000;             % Sampling Frequency in Hz (<-- Set this to match your VI)
F_start_Hz = 1;        % Start Frequency in Hz
F_stop_Hz = 1200;       % Stop Frequency in Hz (Must be < Fs/2)
Desired_RMS = 0.1;     % Desired RMS value
R_realizations = 40;     % Number of realizations (At least 2 required) 

% 1. Convert Hz to Frequency Bins
% Fundamental frequency f0 = Fs / N
f0 = Fs / N;
f_min = round(F_start_Hz / f0);
f_max = round(F_stop_Hz / f0);

% Ensure f_min is at least 1
if f_min < 1; f_min = 1; end

% 2. Generate "Robust Method" Grid (Full Multisine)
% Robust method typically excites ALL integer bins in the band 
Sel = (f_min:f_max)'; 

% 3. Prepare Excitation Matrix (Npp x 2*R)
Sig = zeros(N, R_realizations); 

for r = 1:R_realizations
    % Generate Random Phases (Critical for Robust Method) 
    phases = 2 * pi * rand(size(Sel));
    
    % Construct Frequency Domain Signal
    X = zeros(N, 1);
    X(Sel + 1) = exp(1j * phases); % Map bins to MATLAB indices
    X(N - Sel + 1) = conj(X(Sel + 1)); 
    
    % Inverse FFT to get time domain
    x_time = real(ifft(X));
    
    % --- RMS SCALING ---
    x_time = (x_time / std(x_time)) * Desired_RMS; 

    % Assign to ODD column (Data), EVEN column is 0 (Dummy)
    Sig(:, r ) = x_time;
end

% 4. Save Files
save(fullfile(outputDir, [filename '_Sig_E0_S0.mat']), 'Sig');
save(fullfile(outputDir, [filename '_Sel_E0_S0.mat']), 'Sel');
disp(['Robust signals saved. Freq range: ' num2str(min(Sel)*f0) ' Hz to ' num2str(max(Sel)*f0) ' Hz.']);

% --- PLOTTING ---
figure('Name', 'Robust Method Signal (Hz Scale)');

% Plot Time Domain
subplot(2,1,1);
plot((0:N-1)/Fs, Sig(:,1)); 
xlabel('Time (s)'); ylabel('Amplitude'); 
title(['Time Domain (Realization 1) - RMS: ' num2str(std(Sig(:,1)))]); 
grid on; xlim([0 N/Fs]);

% Plot Frequency Domain
subplot(2,1,2);
F_spectrum = abs(fft(Sig(:,1)));
freq_axis = (0:N-1) * f0; % Frequency axis in Hz
stem(freq_axis, F_spectrum ,'LineWidth', 0.000001, 'Marker', 'none'); 
xlim([0 F_stop_Hz * 1.1]);
xlabel('Frequency (Hz)'); ylabel('Magnitude'); 
title('Frequency Domain (Full Grid Excitation)'); grid on;