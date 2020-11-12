% Micro-Doppler signature of one range bin
close all;

% data
% -------------------------------------------------------------------------
classes = ["walking", "2walking", "running"];          % classes of dataset
for classnum = 3
    class = convertStringsToChars(classes(classnum));
    for setnum = 3
        dataset = [class, num2str(setnum)];       % dataset being run
        
        observedRangeBin = 30;              % observed range bin for micro-Doppler


% read in radar data
% -------------------------------------------------------------------------
folder = fullfile('..\data', class, dataset, 'Radar_Data');   % folder path

load (fullfile(folder, 'CI.mat'))               % load CI file
PRF = CI.PRF;                                   % pulse repetition freq
num_r_bins = CI.num_r_bins;                     % number of range bins
pulses_per_file = CI.pulses_per_file;           % number of pulses per file
fc_MHz = CI.fc_MHz;                             % centre frequency of pulse
start_range_bin = CI.start_range_bin;           % start range bin number
range_res = CI.RR;                              % range resolution

data_files = dir(fullfile(folder, '*.brd'));    % data files
num_data_files = length(data_files);            % number of data files

data = complex(zeros(num_r_bins, ((pulses_per_file+1)*num_data_files)-1, 'single'));   % data matrix

% read first file into data matrix
filepath = fullfile(folder, data_files(1).name);
data(:,1:pulses_per_file) = read_BRDF_binary(filepath, num_r_bins, pulses_per_file);

% read the rest of the files into data matrix
for k = 2: num_data_files
    start_col = k + ((k-1)*pulses_per_file);
    end_col = start_col + pulses_per_file - 1;
    
    filepath = fullfile(folder, data_files(k).name);
    data(:,start_col:end_col) = read_BRDF_binary(filepath, num_r_bins, pulses_per_file);
    data(:,start_col-1) = (data(:,start_col-2)+ data(:,start_col))/2;
end

% stft parameters
% -------------------------------------------------------------------------
fs = PRF;               % sampling frequency

frameLen = 1024;        % frame length
hop = 256;              % hop to next frame
nfft = frameLen;        % number of fft points

% using custom stft function
%--------------------------------------------------------------------------
row = observedRangeBin-start_range_bin+1;   % row in matrix representing range bin
y = data(row, :);                   % data in row

% y = highpass(x, 20, fs);            % filtered data

[s, f, t] = custom_stft(y, frameLen, hop, nfft, fs);    % perform stft

s_dB = 20*log10(abs(s));    % power of s in dB

c = 299792458;              % speed of light (m/s)
lambda = c/(fc_MHz*1e6);    % wavelength of pulse
vel = (f*lambda)/2;         % doppler velocity vector

figure;
imagesc(t,vel,s_dB)
title({['Micro-Doppler of range bin ', num2str(observedRangeBin), ...
    ' (range = ' , num2str(round(observedRangeBin*range_res*100)/100), 'm )'],...
    ['Dataset: ', dataset]})
xlabel('Time (s)')
ylabel('Radial Velocity (m/s)')

hcol = colorbar;
colormap jet
ylabel(hcol, 'Power (dB)')

outputpath = fullfile('micro_doppler_outputs', class, [dataset, '_', num2str(observedRangeBin), '.png']);

set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', outputpath);

ylim([-10 10])
%xlim([70 90])

outputpath = fullfile('micro_doppler_outputs', class, [dataset, '_', num2str(observedRangeBin), '_zoomed.png']);

set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', outputpath);

    end
end
% custom stft function
%--------------------------------------------------------------------------
function [s, f, t] = custom_stft(y, frameLen, hop, nfft, fs)
    y = y(:);                               % transpose to column vector
    
    ylen = length(y);                       % length of vector y
    
    win = hamming(frameLen);                % hamming window
    
    nframes = 1+ fix((ylen-frameLen)/hop);  % number of frames
    
    s = zeros(nfft, nframes);               % initialise stft matrix

    % perform STFT
    for i = 1:nframes
        start = 1+(i-1)*hop;                % starting index of frame
        
        frame = y(start:start+frameLen-1);  % frame
        
        winFrame = frame.*win;              % windowed frame

        % FFT
        fft_winFrame = fftshift(fft(winFrame, nfft));

        % update the stft matrix
        s(:, i) = fft_winFrame;
    end
    
    % calculate time and frequency vectors
    t = (frameLen/2:hop:frameLen/2+(nframes-1)*hop)/fs;
    f = (-nfft/2:1:(nfft/2-1))*fs/nfft;
    
end