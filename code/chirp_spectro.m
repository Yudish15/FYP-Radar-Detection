close all;

% generate chirp
%--------------------------------------------------------------------------
fs = 192000;
start_f = 38000;
stop_f = 40000;
t = 0:1/fs:1;
t1 = 1;

y = chirp(t, start_f, t1, stop_f, 'linear');

figure
plot(t, y)
title('Chirp signal')
xlabel('t')
ylabel('y(t)')

%spectrogram parameters
%--------------------------------------------------------------------------
wlen = 1024;        % window length
num_overlap = 512;  % frame overlap size
hop = wlen - num_overlap;   % hop size
nfft = 1024;        %number of fft points


% using matlab's spectrogram function
%--------------------------------------------------------------------------
figure
spectrogram(y, wlen, num_overlap, nfft, fs, 'yaxis')
colormap jet

% using custom stft function
%--------------------------------------------------------------------------
[stft, f, t] = custom_stft(y, wlen, hop, nfft, fs);
stft_dB = 20*log10(abs(stft));

figure;
imagesc(t,f/1000,stft_dB)
title('Custom Spectrogram')
xlabel('Time(s)')
ylabel('Frequency(kHz)')

hcol = colorbar;
colormap jet
ylabel(hcol, 'Magnitude, dB')


% custom stft function
%--------------------------------------------------------------------------
function [stft, f, t] = custom_stft(y, wlen, hop, nfft, fs)
    y = y(:);   % transpose to column vector
    ylen = length(y);
    
    win = hamming(wlen);
    
    num_frames = 1+ fix((ylen-wlen)/hop);
    
    stft = zeros(nfft, num_frames); % initialise stft matrix
    
    % initialize the indexes
    indx = 0;
    col = 1;

    % perform STFT
    while indx + wlen <= ylen
        % windowing
        yw = y(indx+1:indx+wlen).*win;

        % FFT
        Y = fftshift(fft(yw, nfft));

        % update the stft matrix
        stft(:, col) = Y(1:nfft);

        % update the indexes
        indx = indx + hop;
        col = col + 1;
    end
    
    % calculate time and frequency vectors
    t = (wlen/2:hop:wlen/2+(num_frames-1)*hop)/fs;
    f = (-nfft/2:1:(nfft/2-1))*fs/nfft;
    
end