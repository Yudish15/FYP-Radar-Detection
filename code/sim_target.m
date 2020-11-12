% Target Simulator
close all;

% CFAR parameters
% -------------------------------------------------------------------------
pfa_set = 1e-3;             % probability of false alarm set
N = 16;                     % total number of reference cells
nr = N/2;                   % number of reference cells on one side
ng = 1;                     % number of guard cells on each side

% CA
alpha_ca = (pfa_set.^(-1/N))-1;    % ca cfar constant

% OS
pfa_differ = pfa_set*0.001; % error between pfa_set and pfa_achieved
stat = N*3/4;               % k statistic

alpha_os = 0;               % os cfar constant

% iterative solution to find value of alpha corresponding to pfa set
for alphaVal = 0:0.0001:50      % range of values for alpha
    pfa_comp_os = 1;               % pfa computed from current value of alpha
    for i = 0:stat-1
        pfa_comp_os = pfa_comp_os*(N - i)/(N - i + alphaVal);
    end
    if abs(pfa_comp_os - pfa_set) < pfa_differ
        alpha_os = alphaVal;           
        pfa_differ = abs(pfa_comp_os - pfa_set);
        pfa_achieved_os = pfa_comp_os;    % pfa achieved
    end
end

% SOCA
pfa_differ = pfa_set*0.001; % error between pfa_set and pfa_achieved

alpha_soca = 0;             % soca cfar constant

% iterative solution to find value of alpha corresponding to pfa set
for alphaVal = 0:0.0001:50      % range of values for alpha
    temp = 0;                   % temporary variable used to calculate pfa 
    for k = 0:nr-1
        temp = temp + (factorial(nr-1+k)/(factorial(k)*factorial(nr-1)))*(2+alphaVal)^(-k);
    end
    pfa_comp_soca = 2*((2+alphaVal)^(-nr))*temp;
    if abs(pfa_comp_soca - pfa_set) < pfa_differ
        alpha_soca = alphaVal;           
        pfa_differ = abs(pfa_comp_soca - pfa_set);
        pfa_achieved_soca = pfa_comp_soca;    % pfa achieved
    end
end

% Generate signal with targets and interference
% -------------------------------------------------------------------------
numSamples = 400;       % number of samples

noise = (randn(1,numSamples)+ 1i*randn(1,numSamples))*1/sqrt(2);    % complex gaussian noise

SNR_db_1 = 22;  % snr of 1st target in db
SNR_db_2 = 17;  % snr of 2nd target in db
SNR_db_3 = 19;  % snr of 3rd target in db

% linear snr of targets
snr_1 = 10^(SNR_db_1/20);
snr_2 = 10^(SNR_db_2/20);
snr_3 = 10^(SNR_db_3/20);

signal = noise;
sigLen = length(signal);

% position of targets/interference in signal
pos1 = sigLen/2;
pos2 = (sigLen/2)-4;
pos3 = (sigLen/2)+5;

% insert targets in signal

%signal([pos1]) = signal(pos1)*snr_1;

signal([pos1 pos2 pos3]) = [signal(pos1)*snr_1, signal(pos2)*snr_2, signal(pos3)*snr_3];

signal([pos1 pos2 pos3:end]) = [signal(pos1)*snr_1, signal(pos2)*snr_2, signal(pos3:end)*snr_3];


signal_mag2 = (abs(signal)).^2;  % square law detector

% plot power of signal in dB
figure
plot(10*log10(signal_mag2))
title('CFAR Thresholds')
xlabel('Range bin')
ylabel('Power (dB)')

hold on

% initialise threshold arrays
t_ca = zeros(size(signal_mag2));
t_os = zeros(size(signal_mag2));
t_soca = zeros(size(signal_mag2));

first = 1 + ng + nr;            % first threshold index
last = sigLen - nr - ng;        % last threshold index

% set threshold and count number of false alarms
for i = first:last
    % CA 
    % sum of reference cells
    g_ca = sum(signal_mag2((i-ng-nr):(i-ng-1))) + sum(signal_mag2((i+ng+1):(i+ng+nr)));
    t_ca(i) = g_ca*alpha_ca;     % set CA threshold
    
    % OS
    % sort ref cells in ascending order
    refCells = [(signal_mag2((i-ng-nr):(i-ng-1))), (signal_mag2((i+ng+1):(i+ng+nr)))];
    sortedRefCells = sort(refCells,'ascend');
    g_os = sortedRefCells(stat);        % kth order statistic
    t_os(i) = g_os*alpha_os;            % set OS threshold
    
    % SOCA
    % min(sum lagging window, sum leading window)
    g_soca = min([sum(signal_mag2((i-ng-nr):(i-ng-1))), sum(signal_mag2((i+ng+1):(i+ng+nr)))]);
    t_soca(i) = g_soca*alpha_soca;             % set threshold
    
end

% plot thresholds
plot(10*log10(t_ca))
plot(10*log10(t_os))
plot(10*log10(t_soca))
legend('signal', 'CA-CFAR', 'OS-CFAR', 'SOCA-CFAR')
hold off

ylim([-15 30])
set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', '.\target_outputs\3targets.png');
