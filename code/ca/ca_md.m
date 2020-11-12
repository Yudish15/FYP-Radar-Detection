% CA CFAR on measured data
close all;

% data
% -------------------------------------------------------------------------
classes = ["walking", "2walking", "running"];          % classes of dataset
for classnum = 3
    class = convertStringsToChars(classes(classnum));
    for setnum = 1
        dataset = [class, num2str(setnum)];       % dataset being run
        observedRangeBin = 30;              % observed range bin for micro-Doppler

% read in radar data
% -------------------------------------------------------------------------
folder = fullfile('..\..\data', class, dataset, 'Radar_Data');   % folder path

load (fullfile(folder, 'CI.mat'))               % load CI file
PRF = CI.PRF;                                   % pulse repetition freq
num_r_bins = CI.num_r_bins;                     % number of range bins
pulses_per_file = CI.pulses_per_file;           % number of pulses per file
fc_MHz = CI.fc_MHz;                             % centre frequency of pulse
start_range_bin = CI.start_range_bin;           % start range bin number
range_res = CI.RR;                              % range resolution

data_files = dir(fullfile(folder, '*.brd'));    % data files
num_data_files = length(data_files);            % number of data files
%num_data_files = 1000;                          % number of data files

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

% CA CFAR parameters
% -------------------------------------------------------------------------
pfa_set = 1e-6;         % probability of false alarm set
N = 16;                 % total number of reference cells
nr = N/2;               % number of reference cells on one side
ng = 1;                 % number of guard cells on each side

alpha_ca = (pfa_set.^(-1/N))-1;    % ca cfar constant

% range-Doppler parameters
% -------------------------------------------------------------------------
fs = PRF;               % sampling frequency

frameLen = 2048;        % frame length
hop = 1000;             % hop to next frame
nfft = frameLen;        % number of fft points

% Generate Range-Doppler and perform thresholding
% -------------------------------------------------------------------------
nrow = size(data, 1);                   % number of rows in data
ncol = size(data, 2);                   % number of columns in data

first = 1 + ng + nr;                    % first row to threshold
last = nrow - nr - ng;                  % last row to threshold

nframes = 1 + fix((ncol-frameLen)/hop); % number of frames

win = hamming(frameLen);                % hamming window
wint = win.';                           % transpose win vector to row
winmat = repmat(wint,nrow,1);           % replicate wint row to create a matrix same size as frame

% calculate frequency and range bin vectors
f = (-nfft/2:1:(nfft/2-1)) * fs/nfft;                       % frequency range
r = (start_range_bin:start_range_bin + num_r_bins - 1)*range_res;     % range

t = (frameLen/2:hop:frameLen/2+(nframes-1)*hop)/fs;         % time vector

detMat = NaN(nframes, num_r_bins);              % detection matrix for whole data

observedRow = observedRangeBin-start_range_bin+1;   % row in matrix representing observed range bin

observedRangeBin_fft_mat = zeros(nframes, nfft);    % fft matrix for observed range bin
observedRangeBin_det_mat = zeros(nframes, nfft);    % detection matrix for observed range bin

for i = 1:nframes
    
    disp = (i == 100) || (i == 500) || (i == 1000);
    
    start = 1+(i-1)*hop;                        % start column index of frame

    frame = data(:, start : start+frameLen-1);  % frame

    winframe = frame.*winmat;                   % windowed frame

    fft_winframe = fftshift(fft(winframe, [], 2), 2);   % fft of windowed frame
    
    fft_winframe_mag2 = (abs(fft_winframe)).^2;         % fft power
    
    observedRangeBin_fft_mat(i, :) = fft_winframe_mag2(observedRow,:);  
    
%     fft_winframe_mag2(:,(frameLen/2)-19:(frameLen/2)+21) = 0;

    if disp
        % plot fft power
        figure;
        imagesc(f/1000, r, 10*log10(fft_winframe_mag2))
        title({['CA-CFAR Range-Doppler map (frame: ', num2str(i), ')'],['Dataset: ', dataset]})
        xlabel('Doppler frequency (kHz)')
        ylabel('Range (m)')

        hcol = colorbar;
        colormap jet
        ylabel(hcol, 'Power (dB)')
        hold on;
    end
   
    t_ca = zeros(size(fft_winframe_mag2));   % initialise threshold array for current frame
    
    % set threshold and count number of detections
    for j = first:last                      % 'first' to 'last' row
        numDet = 0;                         % number of detections for range bin in current frame
        
        for k = 1:length(fft_winframe_mag2(j,:))    % columns
            % sum reference cells
            g_ca = sum(fft_winframe_mag2((j-ng-nr):(j-ng-1),k)) + sum(fft_winframe_mag2((j+ng+1):(j+ng+nr),k));
            % set threshold for cell under test
            t_ca(j, k) = g_ca*alpha_ca;

            % check for detection
            if t_ca(j,k) < fft_winframe_mag2(j,k)
                numDet = numDet+1;
                if disp
                    plot(f(k)/1000,r(j),'kx', 'MarkerSize',6, 'LineWidth',1.2);
                end
            end
        end
        if numDet > 2
            detMat(i,j) = r(j);
        end
    end
    
    observedRangeBin_det_mat(i, :) = t_ca(observedRow,:) < fft_winframe_mag2(observedRow,:);
    
    if disp
        hold off
        outputpath = fullfile('outputs', class, dataset, ['range_doppler_', num2str(i), '.png']);
        set(gcf,'PaperPosition',[0 0 16 10])
        print(gcf, '-dpng', outputpath);
        
        xlim([-1 1])
        outputpath = fullfile('outputs', class, dataset, ['range_doppler_', num2str(i), '_zoomed.png']);
        set(gcf,'PaperPosition',[0 0 16 10])
        print(gcf, '-dpng', outputpath);
    end

end

% plot target range over time
figure
plot(t, detMat, 'k')          % plots each column of the matrix as a series
ylim([r(1) r(end)])

title({'CA-CFAR Target Range vs Time', ['(Dataset: ', dataset, ')']})
xlabel('Time (s)')
ylabel('Range (m)')

outputpath = fullfile('outputs', class, dataset, ['location_vs_time.png']);
set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', outputpath);

% plot spectrogram of observed range bin with detections
figure
imagesc(t, f/1000, 10*log10(observedRangeBin_fft_mat.'))
title({['CA-CFAR Detections on Range Bin ', num2str(observedRangeBin),...
    ' (range = ', num2str(round(r(observedRow)*100)/100), 'm )'],['Dataset: ', dataset]})
xlabel('Time (s)')
ylabel('Doppler frequency (kHz)')

hcol = colorbar;
colormap jet
ylabel(hcol, 'Power (dB)')
hold on;

[rpos, cpos] = find (observedRangeBin_det_mat.');
plot(t(cpos),f(rpos)/1000, 'kx', 'MarkerSize',6, 'LineWidth',1.2);

hold off

outputpath = fullfile('outputs', class, dataset, ...
    ['rb', num2str(observedRangeBin), '_n', num2str(N), '_ng', num2str(ng), '_', num2str(pfa_set),'.png']);
set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', outputpath);

    end
end
