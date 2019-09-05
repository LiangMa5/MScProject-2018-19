% RSMA partial CSIT SR vs SNR
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% MISO: Nr=1 / Nt = N_user =3
% Rate region for SDMA, RSMA(1 layer)
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc;clear all;

Nr = 1; N_user = 2; Nt = 2;
SNRdB = 5:5:30;  %SNR in dB
SNR = 10.^(SNRdB/10);
weight = ones(1,N_user);

%accuracy of convergence
tolerance = 1e-6;
M = 100;

c = 1; d = 0.3; %error channel power = c*SNR^(-d)
power_err = c*SNR.^(-d);

Rth = [0.15; 0.2; 0.25; 0.3; 0.3; 0.3];
rng(1)
%% rate region simulation
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));
% H_err generation
for i_m = 1:M
    H_err(:,:,i_m) = sqrt(0.5)*randn(N_user,Nt)+1i*sqrt(0.5)*randn(N_user,Nt);
end
%%
ite = 100;
%H_est generation
for i0 = 1:ite
    H_est(:,:,i0) = sqrt(0.5)*randn(N_user,Nt)+1i*sqrt(0.5)*randn(N_user,Nt); % random channel
end

% 100 realizations
for i1 = 1:ite
    tic
    for i_snr = 1:length(SNRdB)
        H = sqrt(1-power_err(i_snr))*H_est(:,:,i1) + sqrt(power_err(i_snr))*H_err;
        Rate_SDMA(i_snr,i1) = SDMA_Rate(H_est(:,:,i1),H,SNRdB(i_snr),weight,tolerance,M,c,d,Nt,N_user,Rth(i_snr));
        Rate_RS1layer(i_snr,i1) = RS_1layer_Rate(H_est(:,:,i1),H,SNRdB(i_snr),weight,tolerance,M,c,d,Nt,N_user,Rth(i_snr));
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

