% RSMA CoMP JT SR vs SNR
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% SISO: Nt=1, Nr=1 / 3 base stations / 3 user
% WSR vs SNR for SDMA, RSMA(1 layer)
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc; clear all;

% channel
Nr = 1; N_bs = 3;%number of base station
N_user = 3;
SNRdB = 5:5:30;  %SNR in dB
Rth = [0:12; 0:18; 0:24; 0:3; 0:3; 0:3];
%user weights
weight = [1,1,1];


%accuracy of convergence
tolerance = 1e-6;

alpha = 1;
beta = 1;
ab = [1,alpha,0;alpha*beta,beta,alpha*beta;0,alpha,1]; %entry: UE*BS


%% WSR
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));


for i1 = 1:100
    tic
    Hran = sqrt(0.5)*randn(Nr,N_bs,N_user)+1i*sqrt(0.5)*randn(Nr,N_bs,N_user); % random channels
    for i0 = 1:N_user
        H(:,:,i0) = sqrt(ab(i0,:)).* Hran(:,:,i0); % apply channel strength
    end
    
    parfor i_snr = 1:length(SNRdB)
        Rate_SDMA(i_snr,i1) = SDMA_Rate(H,SNRdB(i_snr),weight,tolerance,Rth(i_snr));
        Rate_RS1layer(i_snr,i1) = RS_1layer_Rate(H,SNRdB(i_snr),weight,tolerance,Rth(i_snr));
    end
    
    save('Rate_SDMA.mat','Rate_SDMA');
    save('Rate_RS1layer.mat','Rate_RS1layer');
    fprintf('loop %d done   ',i1);
    toc
end

%%
figure (1)
plot(SNRdB,mean(Rate_SDMA,2),'-.','LineWidth',2.5); hold on;grid on
plot(SNRdB,mean(Rate_RS1layer,2),'o-','LineWidth',2.5); grid on
xlabel('SNR (dB)');
ylabel('WSR (bits/s/Hz)');
legend('SDMA','RSMA');

