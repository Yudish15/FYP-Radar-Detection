% 1D CA-CFAR Simulation
close all;

% CA CFAR parameters
% -------------------------------------------------------------------------
pfa_set = 1e-2;         % probability of false alarm set
N = 16;                 % total number of reference cells
nr = N/2;               % number of reference cells on one side
ng = 1;                 % number of guard cells on each side

alpha_ca = (pfa_set.^(-1/N))-1;    % ca cfar constant

% Generate complex gaussian noise
% -------------------------------------------------------------------------
numSamples = 2e5;       % number of samples

noise = (randn(1,numSamples)+ 1i*randn(1,numSamples))*1/sqrt(2);    % complex gaussian noise

%plot(noise, '.')

noise_mag2 = (abs(noise)).^2;       % magnitude of noise squared (square law detector)

% plot magnitude squared of noise
figure
plot(noise_mag2)
title('1D CA CFAR Thresholding')
xlabel('Samples')
ylabel('Power')

hold on

t_ca = zeros(size(noise_mag2));         % initialise threshold array
first = 1 + ng + nr;                    % first threshold index
last = length(noise) - nr - ng;         % last threshold index

numFA = 0;                              % number of false alarms

% set threshold and count number of false alarms
for i = first:last
    g_ca = sum(noise_mag2((i-ng-nr):(i-ng-1))) + sum(noise_mag2((i+ng+1):(i+ng+nr)));    % sum of reference cells
    t_ca(i) = g_ca*alpha_ca;     % set threshold
    
    % check if false alarm
    if t_ca(i) < noise_mag2(i)
        numFA = numFA + 1;
    end
    
end

% plot threshold
plot (t_ca)
legend('noise', 'threshold')
hold off

set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', '.\outputs\ca_sim1_full.png');

xlim([1000 1500])
set(gcf,'PaperPosition',[0 0 16 10])
print(gcf, '-dpng', '.\outputs\ca_sim1.png');

pfa_obtained = numFA/numSamples              % probability of false alarm obtained
pfa_error = (abs(pfa_set-pfa_obtained)/pfa_set)*100     % error in pfa

fid = fopen('.\outputs\ca_sim1.txt','a');
fprintf(fid, 'PFA set: %.4d \nPFA obtained: %.4d \nPFA error: %.2f%%\n', pfa_set, pfa_obtained, pfa_error);
fclose(fid);