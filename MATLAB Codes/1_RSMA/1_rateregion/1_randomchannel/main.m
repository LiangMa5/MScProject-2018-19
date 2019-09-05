% RSMA Rate region (random channel)
% Implemented algorithm in the programme is adopted from the paper:
% Rate-splitting multiple access for downlink communication 
% systems: bridging, generalizing, and outperforming SDMA and NOMA.
% Equations can be found in the above paper using the indicated number.

% MISO: Nt=4, Nr=1 / two user
% Rate region for SDMA, NOMA and RSMA
% aim: find the best precoder which maximizes the rate

%% parameter setting
clc; clear all;

% channel
Nr = 1; N_t = 2;
N_user = 2;
SNRdB = 20;  %SNR in dB

%user weights
u1 = 1;
u2 = 10.^[-3 -1:0.05:1 3];

%accuracy of convergence
tolerance = 1e-6;

sigma2 = 1;
%% 100 random channel realizations
clk = fix(clock); fprintf('Start time is %d:%d  \n', clk(4),clk(5));
for i1 = 1:100
    tic
    H = sqrt(0.5)*randn(Nr,N_t,N_user)+1i*sqrt(0.5)*randn(Nr,N_t,N_user); % random channels
    H(:,:,2) = sqrt(sigma2)*H(:,:,2);
    parfor i_u2 = 1:length(u2)
        weight = [u1,u2(i_u2)];
        Rate_SDMA(i_u2,:) = SDMA_Rate(H,SNRdB,weight,tolerance);
        [Rate_order1n(i_u2,:),Rate_order2n(i_u2,:)] = NOMA_Rate(H,SNRdB,weight,100*tolerance);
        [Rate_order1(i_u2,:),Rate_order2(i_u2,:)] = RS_Rate(H,SNRdB,weight,tolerance);
    end
    
    Rate_SDMA1(:,i1) = Rate_SDMA(:,1);
    Rate_SDMA2(:,i1) = Rate_SDMA(:,2);
    
    Rate_NOMA1(:,i1) = [Rate_order1n(:,1);Rate_order2n(:,1)];
    Rate_NOMA2(:,i1) = [Rate_order1n(:,2);Rate_order2n(:,2)];
    
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
xx = floor(10*x1)/10;
indexmin = find(xx==0);
[~,indexmax] = max(x1);

figure (1)
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),'-.','LineWidth',3); hold on;grid on
xlabel('{\it{R_{total,1}}} (bits/s/Hz)');
ylabel('{\it{R_{total,2}}} (bits/s/Hz)');


%% NOMA plot
v = mean(Rate_NOMA1,2);
w = mean(Rate_NOMA2,2);
save('v.mat','v');save('w.mat','w');

k = convhull(v,w);
x1 = v(k);
y1 = w(k);
xx = floor(x1);
indexmin = find(xx==0);
[~,indexmax] = max(x1);
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),':','LineWidth',2.5); hold on;



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
legend('SDMA','NOMA','RSMA');

