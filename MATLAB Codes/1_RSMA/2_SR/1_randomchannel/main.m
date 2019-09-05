% RSMA SNR SR vs SNR
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% MISO: Nt=3, Nr=1 / 3 user
% WSR vs SNR for SDMA, NOMA and RSMA
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc; clear all;


Nr = 1; N_t = 3;%number of base station
N_user = 3;

SNRdB = 0:5:30;  %SNR in dB
Rth = [0.1; 0.15; 0.2; 0.25; 0.3; 0.3; 0.3];
%user weights
weight = [1,1,1];
%accuracy of convergence
tolerance = 1e-6;


%% WSR
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));

for i1 = 1:100
    tic
    H = sqrt(0.5)*randn(Nr,N_t,N_user)+1i*sqrt(0.5)*randn(Nr,N_t,N_user); % random channels
    for i_snr = 1:length(SNRdB)
        Rate_SDMA(i_snr,i1) = SDMA_Rate(H,SNRdB(i_snr),weight,tolerance,Rth(i_snr));
        Rate_NOMA(i_snr,i1)  = NOMA_Rate(H,SNRdB(i_snr),weight,tolerance,Rth(i_snr));
        Rate_RS(i_snr,i1)  = RS_Rate(H,SNRdB(i_snr),weight,tolerance,Rth(i_snr));
        Rate_RS1layer(i_snr,i1) = RS_1layer_Rate(H,SNRdB(i_snr),weight,tolerance,Rth(i_snr));
    end
    
    save('Rate_SDMA.mat','Rate_SDMA');
    save('Rate_NOMA.mat','Rate_NOMA');
    save('Rate_RS.mat','Rate_RS');
    save('Rate_RS1layer.mat','Rate_RS1layer');
    
    fprintf('loop %d done   ',i1);
    toc
end

%%
figure (1)
plot(SNRdB,mean(Rate_SDMA,2),'-.','LineWidth',2.5); hold on;
plot(SNRdB,mean(Rate_NOMA,2),':','LineWidth',2.5); hold on;
plot(SNRdB,mean(Rate_RS,2),'*-','LineWidth',2.5); hold on;
plot(SNRdB,mean(Rate_RS1layer,2),'o-','LineWidth',2.5); grid on
xlabel('SNR (dB)');
ylabel('WSR (bits/s/Hz)');
legend('SDMA','NOMA','RS','RS 1-layer');

