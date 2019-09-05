% RSMA CoMP JT partial CSIT channel disparity
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% SISO: Nt=1, Nr=1 / 3 user / 3 BS
% Rate region for SDMA, RSMA (1 layer)
% aim: find the best precoder which maximizes the rate

%% parameter setting
clear all;

Nr = 1; N_user = 3; N_bs = 3; % Nt appears in this folder is considered as N_bs
SNRdB = 20;  %SNR in dB
SNR = 10.^(SNRdB/10);
weight = ones(1,N_user);

%accuracy of convergence  
tolerance = 1e-6;
M = 100;

c = 1; d = 0.6; %error channel power = c*SNR^(-d)
power_err = c*SNR.^(-d);

Rth = 0.1;

alpha = 0.1:0.1:1;  beta = 0.1;
for i_a = 1:length(alpha)
    ab(:,:,i_a) = [1,alpha(i_a),0;alpha(i_a)*beta,beta,alpha(i_a)*beta;0,alpha(i_a),1]; %entry: UE*BS
end

%% rate region simulation
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));
% H_err generation
for i_m = 1:M
    H_err(:,:,i_m) = sqrt(0.5)*randn(N_user,N_bs)+1i*sqrt(0.5)*randn(N_user,N_bs);
end

ite = 3;
%H_est generation
for i_0 = 1:ite
    H_est0(:,:,i_0) = sqrt(0.5)*randn(N_user,N_bs)+1i*sqrt(0.5)*randn(N_user,N_bs); % random channel 
end

% 100 realizations
for i1 = 1:ite
    tic
    parfor i_ab = 1:length(alpha)
        H_est = sqrt(ab(:,:,i_ab)).*H_est0(:,:,i1);     
        H = sqrt(1-power_err)*H_est + sqrt(power_err)*H_err;
        Rate_SDMA(i_ab,i1) = SDMA_Rate(H_est,H,SNRdB,weight,tolerance,M,c,d,N_bs,N_user,Rth);
        Rate_RS1layer(i_ab,i1) = RS_1layer_Rate(H_est,H,SNRdB,weight,tolerance,M,c,d,N_bs,N_user,Rth);
    end
   
    save('Rate_SDMA.mat','Rate_SDMA');
    save('Rate_RS1layer.mat','Rate_RS1layer');
    fprintf('loop %d done   ',i1);
    toc
end



%%
figure (1)
plot(alpha,mean(Rate_SDMA,2),'-.','LineWidth',2.5); hold on;grid on
plot(alpha,mean(Rate_RS1layer,2),'o-','LineWidth',2.5); grid on
xlabel('SNR (dB)');
ylabel('WSR (bits/s/Hz)');
legend('SDMA','RSMA');

