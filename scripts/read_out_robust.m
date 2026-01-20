clear; clc;
%% ------- Output display ----------
Fs = 4000;


[u, y, sel, sig, ~, ~] = acquisition("results_4_fast/", "outOddM");

  % remove the second realization (useless)
    u = u(:, 1);
    y = y(:, 1);
    N = length(u) /10; % Calculate the number of samples per period
  % transient removal
    assert(length(u) == 10*N);
    u = u(N+1:end);
    y = y(N+1:end);

   
U = fft(u);    
Y = fft(y);

figure;
subplot(211);
plot((0:10*N-N-1)*Fs/N, db(U), '.', LineWidth=0.5);
%xlim([1/fs maxFreq]);
ylim([min(db(U)) 1.5*max(db(U))]);
title('Input spectrum');
xlabel('Frequency [Hz]');
ylabel('Amplitude [dBV]');

subplot(212);
plot((0:10*N-N-1)*Fs/N, db(Y), '.', LineWidth=0.1);
%xlim([1/fs maxFreq]);
title('Output spectrum');
xlabel('Frequency [Hz]');
ylabel('Amplitude [dBV]');
