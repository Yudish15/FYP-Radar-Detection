% 2D SOCA-CFAR Simulation
close all;

% SOCA CFAR parameters
% -------------------------------------------------------------------------
pfa_set = 1e-3;             % probability of false alarm set
pfa_differ = pfa_set*0.001; % error between pfa_set and pfa_achieved
N = 16;                     % total number of reference cells
nr = N/2;                   % number of reference cells on one side
ng = 1;                     % number of guard cells on each side

alpha_soca = 0;                  % soca cfar constant

% iterative solution to find value of alpha corresponding to pfa set
for alphaVal = 0:0.0001:50      % range of values for alpha
    temp = 0;                   % temporary variable used to calculate pfa 
    for k = 0:nr-1
        temp = temp + (factorial(nr-1+k)/(factorial(k)*factorial(nr-1)))*(2+alphaVal)^(-k);
    end
    pfa_comp = 2*((2+alphaVal)^(-nr))*temp;
    if abs(pfa_comp - pfa_set) < pfa_differ
        alpha_soca = alphaVal;           
        pfa_differ = abs(pfa_comp - pfa_set);
        pfa_achieved = pfa_comp;    % pfa achieved
    end
end

% Generate complex gaussian noise matrix
% -------------------------------------------------------------------------
numRangeBins = 51;      % number of range bins
numCol = 20000;        % number of columns

noiseMatrix = (randn(numRangeBins,numCol)+1i*randn(numRangeBins,numCol))*1/sqrt(2); % complex noise matrix

data = noiseMatrix;

% range-Doppler parameters
% -------------------------------------------------------------------------
fs = 10000;             % sampling frequency

frameLen = 1024;        % frame length
hop = 1024;             % hop to next frame
nfft = frameLen;        % number of fft points

% Generate Range-Doppler and perform thresholding
% -------------------------------------------------------------------------
nrow = size(data, 1);                   % number of rows in data
ncol = size(data, 2);                   % number of columns in data

first = 1 + ng + nr;                    % first row to threshold
last = nrow - nr - ng;                  % last row to threshold
numFA = 0;                              % number of false alarms

nframes = 1 + fix((ncol-frameLen)/hop); % number of frames

win = hamming(frameLen);                % hamming window
wint = win.';                           % transpose win vector to row
winmat = repmat(wint,nrow,1);           % replicate wint row to create a matrix same size as frame

% calculate frequency and range bin vectors
f = (-nfft/2:1:(nfft/2-1)) * fs/nfft;   % frequency range
r = 1:nrow;                             % range bins

fig1 = figure();

for i = 1:nframes
    start = 1+(i-1)*hop;                        % start column index of frame

    frame = data(:, start : start+frameLen-1);  % frame

    winframe = frame.*winmat;                   % windowed frame

    fft_winframe = fftshift(fft(winframe, [], 2), 2);   % fft of windowed frame
    
    fft_winframe_mag2 = (abs(fft_winframe)).^2;         % fft power
    
    % plot fft power
    figure(fig1)
    imagesc(f/1000, r, 10*log10(fft_winframe_mag2))
    title(['Frame ' num2str(i)])
    xlabel('Doppler frequency(kHz)')
    ylabel('Range bin')

    hcol = colorbar;
    colormap jet
    ylabel(hcol, 'Power (dB)')
    hold on;
   
    t_soca = zeros(size(fft_winframe_mag2));   % initialise threshold array for current frame
    
    % set threshold and count number of false alarms
    for j = first:last                              % 'first' to 'last' row
        for k = 1:length(fft_winframe_mag2(j,:))    % columns
            % g = min(sum lagging window, sum leading window)
            g_soca = min([sum(fft_winframe_mag2((j-ng-nr):(j-ng-1),k)), sum(fft_winframe_mag2((j+ng+1):(j+ng+nr),k))]);
            % set threshold for cell under test
            t_soca(j, k) = g_soca*alpha_soca;

            % check for false alarms
            if t_soca(j,k) < fft_winframe_mag2(j,k)
                numFA = numFA + 1;
                plot(f(k)/1000,r(j),'kx', 'MarkerSize',8, 'LineWidth',1.5);
            end
        end
    end
    
    hold off
    drawnow
end

set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', '.\outputs\soca_sim2.png');

numRows = last-first+1;                 % number of rows for which threshold was calculated
numSamples = numRows*numCol;            % number of samples for which threshold was calculated

pfa_obtained = numFA/numSamples              % probability of false alarm obtained
pfa_error = (abs(pfa_set-pfa_obtained)/pfa_set)*100     % error in pfa

fid = fopen('.\outputs\soca_sim2.txt','a');
fprintf(fid, 'PFA set: %.4d \nPFA obtained: %.4d \nPFA error: %.2f%%\n', pfa_set, pfa_obtained, pfa_error);
fclose(fid);