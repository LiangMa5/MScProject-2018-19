% RSMA partial CSIT rate region
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% MISO: Nt=2, Nr=1 / two user
% Rate region for SDMA, RSMA
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc; clear all;

Nr = 1; Nt = 2; N_user = 2;
SNRdB = 20;  %SNR in dB
SNR = 10.^(SNRdB/10);
u1 = 1;
u2 = 10.^[-3 -1:0.05:1 3];

%accuracy of convergence
tolerance = 1e-6;
M = 100;

c = 1; d = 0.6; %error channel power = c*SNR^(-d)
power_err = c*SNR.^(-d);
%% rate region simulation
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));
% H_err generation
for i_m = 1:M
    H_err(:,:,i_m) = sqrt(0.5)*randn(N_user,Nt)+1i*sqrt(0.5)*randn(N_user,Nt);
end

ite = 100;
%H_est generation
for i0 = 1:ite
    H_est(:,:,i0) = sqrt(0.5)*randn(N_user,Nt)+1i*sqrt(0.5)*randn(N_user,Nt); % random channel
end

% 100 realizations
for i1 = 1:ite
    tic
    H = sqrt(1-power_err)*H_est(:,:,i1) + sqrt(power_err)*H_err;
    for i_u2 = 1:length(u2)
        weight = [u1,u2(i_u2)];
        Rate_SDMA(i_u2,:) = SDMA_Rate(H_est(:,:,i1),H,SNRdB,weight,tolerance,M,c,d,Nt,N_user);
        [Rate_order1(i_u2,:),Rate_order2(i_u2,:)] = RS_Rate(H_est(:,:,i1),H,SNRdB,weight,tolerance,M,c,d,Nt,N_user);
    end

    Rate_SDMA1(:,i1) = Rate_SDMA(:,1);
    Rate_SDMA2(:,i1) = Rate_SDMA(:,2);
    
    Rate_RS1(:,i1) = [Rate_order1(:,1);Rate_order2(:,1)];
    Rate_RS2(:,i1) = [Rate_order1(:,2);Rate_order2(:,2)];
    
    fprintf('loop %d done   ',i1);
    toc
end


%% SDMA plot
x = mean(Rate_SDMA1,2); %rate of user 1
y = mean(Rate_SDMA2,2); %rate of user 2
save('x.mat','x');save('y.mat','y');
k = convhull(x,y);
x1 = x(k);
y1 = y(k);
xx = floor(x1);
indexmin = find(xx==0);
[~,indexmax] = max(x1);

figure (1)
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),'-.','LineWidth',3); hold on;grid on
xlabel('{\it{R_{total,1}}} (bits/s/Hz)');
ylabel('{\it{R_{total,2}}} (bits/s/Hz)');


%% RS plot
t = mean(Rate_RS1,2);
z = mean(Rate_RS2,2);
save('t.mat','t');save('z.mat','z');

k = convhull(t,z);
x1 = t(k);
y1 = z(k);
xx = floor(x1);
indexmin = find(xx==0);
[~,indexmax] = max(x1);
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),'*-','LineWidth',3);
legend('SDMA','RSMA');


