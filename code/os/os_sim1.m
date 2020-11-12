% 1D OS-CFAR Simulation
close all;

% OS CFAR parameters
% -------------------------------------------------------------------------
pfa_set = 1e-3;             % probability of false alarm set
pfa_differ = pfa_set*0.001; % error between pfa_set and pfa_achieved
N = 16;                     % total number of reference cells
stat = N*3/4;               % order statistic
nr = N/2;                   % number of reference cells on one side
ng = 1;                     % number of guard cells on each side

alpha_os = 0;                  % os cfar constant

% iterative solution to find value of alpha corresponding to pfa set
for alphaVal = 0:0.0001:50      % range of values for alpha
    pfa_comp = 1;               % pfa computed from current value of alpha
    for i = 0:stat-1
        pfa_comp = pfa_comp*(N - i)/(N - i + alphaVal);
    end
    if abs(pfa_comp - pfa_set) < pfa_differ
        alpha_os = alphaVal;           
        pfa_differ = abs(pfa_comp - pfa_set);
        pfa_achieved = pfa_comp;    % pfa achieved
    end
end

% Generate complex gaussian noise
% -------------------------------------------------------------------------
numSamples = 2e5;           % number of samples

noise = (randn(1,numSamples)+ 1i*randn(1,numSamples))*1/sqrt(2);    % complex gaussian noise

%plot(noise, '.')

noise_mag2 = (abs(noise)).^2;       % magnitude of noise squared (square law detector)

% plot magnitude squared of noise
figure
plot(noise_mag2)
title('1D OS CFAR Thresholding')
xlabel('Samples')
ylabel('Power')

hold on

t_os = zeros(size(noise_mag2));          % initialise threshold array
first = 1 + ng + nr;                    % first threshold index
last = length(noise) - nr - ng;         % last threshold index

numFA = 0;                              % number of false alarms

% set threshold and count number of false alarms
for i = first:last
    refCells = [(noise_mag2((i-ng-nr):(i-ng-1))), (noise_mag2((i+ng+1):(i+ng+nr)))];    % reference cells
    sortedRefCells = sort(refCells,'ascend');   % sort refCells in ascending order
    g_os = sortedRefCells(stat);        % kth order statistic
    t_os(i) = g_os*alpha_os;             % set threshold
    
    % check if false alarm
    if t_os(i) < noise_mag2(i)
        numFA = numFA + 1;
    end
    
end

% plot threshold
plot (t_os)
legend('noise', 'threshold')
hold off

set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', '.\outputs\os_sim1_full.png');

xlim([1000 1500])
set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', '.\outputs\os_sim1.png');

pfa_obtained = numFA/numSamples              % probability of false alarm obtained
pfa_error = (abs(pfa_set-pfa_obtained)/pfa_set)*100         % error in pfa

fid = fopen('.\outputs\os_sim1.txt','a');
fprintf(fid, 'PFA set: %.4d \nPFA obtained: %.4d \nPFA error: %.2f%%\n', pfa_set, pfa_obtained, pfa_error);
fclose(fid);