% RSMA CoMP JT partial CSIT SR vs SNR
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% SISO: Nt=1, Nr=1 / 3 user / 3 BS
% Rate region for SDMA, RSMA (1 layer)
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc;clear all;

Nr = 1; N_user = 2; N_bs = 2; % Nt appears in this folder is considered as N_bs
SNRdB = 20:5:30;  %SNR in dB
SNR = 10.^(SNRdB/10);
weight = ones(1,N_user);

%accuracy of convergence  
tolerance = 1e-6;
M = 100;

alpha = 0.5;
beta = 0.5;
%ab = [1,alpha,0;alpha*beta,beta,alpha*beta;0,alpha,1]; %entry: UE*BS
%ab = [1,alpha,0,0,0;alpha*beta,beta,alpha*beta,0,0;0,alpha,1,alpha,0;0,0,alpha*beta,beta,alpha*beta;0,0,0,alpha,1];
ab = [1,alpha;alpha*beta,beta]; %entry: UE*BS
c = 1; d = 0.6; %error channel power = c*SNR^(-d)
power_err = c*SNR.^(-d);

Rth = [0.04; 0.06; 0.08; 0.1; 0.1; 0.1]; 
% rate region simulation
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));
% H_err generation
for i_m = 1:M
    H_err(:,:,i_m) = sqrt(0.5)*randn(N_user,N_bs)+1i*sqrt(0.5)*randn(N_user,N_bs);
end
%%
ite = 100;
%H_est generation
for i0 = 1:ite
    H_est0(:,:,i0) = sqrt(0.5)*randn(N_user,N_bs)+1i*sqrt(0.5)*randn(N_user,N_bs); % random channel
    H_est(:,:,i0) = sqrt(ab).*H_est0(:,:,i0);
end

% 100 realizations
for i1 = 1:ite
    tic
    for i_snr = 1:length(SNRdB)
        H = sqrt(1-power_err(i_snr))*H_est(:,:,i1) + sqrt(power_err(i_snr))*H_err;
        Rate_SDMA(i_snr,i1) = SDMA_Rate(H_est(:,:,i1),H,SNRdB(i_snr),weight,tolerance,M,c,d,N_bs,N_user,Rth(i_snr))
        Rate_RS1layer(i_snr,i1) = RS_1layer_Rate(H_est(:,:,i1),H,SNRdB(i_snr),weight,tolerance,M,c,d,N_bs,N_user,Rth(i_snr));
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

