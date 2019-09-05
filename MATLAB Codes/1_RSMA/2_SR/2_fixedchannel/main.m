% RSMA SR vs SNR
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% MISO: Nt=4, Nr=1 / 3 user
% WSR vs SNR for SDMA, NOMA and RSMA
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc; clear all;

% channel
theta2 = pi/9; %channel angle for h2
theta3 = 2*theta2;
gamma1 = 1; %channel bias
gamma2 = 0.3;
H(:,:,1) = [1,1,1,1];  %channel 1
H(:,:,2) = gamma1*[1,exp(1i*theta2), exp(1i*2*theta2), exp(1i*3*theta2)];
H(:,:,3) = gamma2*[1,exp(1i*theta3), exp(1i*2*theta3), exp(1i*3*theta3)];


% H(:,:,1) = [1,1];  %channel 1
% H(:,:,2) = gamma1*[1,exp(1i*theta2)];
% H(:,:,3) = gamma2*[1,exp(1i*theta3)];


SNRdB = 0:5:30;  %SNR in dB
Rth = [0:02; 0:04; 0:06; 0:08; 0:1; 0:1; 0:1];
%user weights
weight = [0.2,0.3,0.5];
%accuracy of convergence
tolerance = 1e-6;


%% WSR
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));
for i1 = 1:length(SNRdB)
    tic
    Rate_SDMA(i1,:) = SDMA_Rate(H,SNRdB(i1),weight,tolerance,Rth(i1));
    Rate_NOMA(i1,:)  = NOMA_Rate(H,SNRdB(i1),weight,tolerance,Rth(i1));
    Rate_RS(i1,:)  = RS_Rate(H,SNRdB(i1),weight,tolerance,Rth(i1));
    Rate_RS1layer(i1,:) = RS_1layer_Rate(H,SNRdB(i1),weight,tolerance,Rth(i1));
    
    fprintf('loop %d done   ',i1);
    toc
end

save('Rate_SDMA.mat','Rate_SDMA');
save('Rate_NOMA.mat','Rate_NOMA');
save('Rate_RS.mat','Rate_RS');
save('Rate_RS1layer.mat','Rate_RS1layer');

%%
figure (1)
plot(SNRdB,Rate_SDMA,'-.','LineWidth',2.5); hold on;
plot(SNRdB,Rate_NOMA,':','LineWidth',2.5); hold on;
plot(SNRdB,Rate_RS,'*-','LineWidth',2.5); hold on;
plot(SNRdB,Rate_RS1layer,'o-','LineWidth',2.5); grid on
xlabel('SNR (dB)');
ylabel('WSR (bits/s/Hz)');
legend('SDMA','NOMA','RS','RS 1-layer');

