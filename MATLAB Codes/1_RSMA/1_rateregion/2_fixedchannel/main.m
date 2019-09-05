% RSMA Rate region
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
H(:,:,1) = [1,1,1,1];  %channel 1
%H(:,:,1) = [1,1];
theta = pi/9; %channel angle for h2
gamma = 1; %channel bias
H(:,:,2) = gamma*[1,exp(1i*theta), exp(1i*2*theta), exp(1i*3*theta)];
%H(:,:,2) = gamma*[1,exp(1i*theta)];

SNRdB = 20;  %SNR in dB

%user weights
u1 = 1;
u2 = 10.^[-3 -1:0.05:1 3];

%accuracy of convergence
tolerance = 1e-6;

%% rate region simulation
clk = fix(clock);  fprintf('Start time is %d:%d  \n', clk(4),clk(5));

parfor i_u2 = 1:length(u2)
    tic
    weight = [u1,u2(i_u2)];
    Rate_SDMA(i_u2,:) = SDMA_Rate(H,SNRdB,weight,tolerance);
    [Rate_order1n(i_u2,:),Rate_order2n(i_u2,:)] = NOMA_Rate(H,SNRdB,weight,tolerance);
    [Rate_order1(i_u2,:),Rate_order2(i_u2,:)] = RS_Rate(H,SNRdB,weight,tolerance);
    fprintf('loop %d done   ',i_u2);
    toc
end


%% SDMA plot
x = Rate_SDMA(:,1); %rate of user 1
y = Rate_SDMA(:,2); %rate of user 2
save('x.mat','x');save('y.mat','y');

k = convhull(x,y);
x1 = x(k);
y1 = y(k);
xx = floor(x1);
indexmin = find(xx==0);
[~,indexmax] = max(x1);
figure (1)
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),'-.','LineWidth',2.5); hold on;grid on

%% NOMA plot
v = [Rate_order1n(:,1);Rate_order2n(:,1)];
w = [Rate_order1n(:,2);Rate_order2n(:,2)];
save('v.mat','v');save('w.mat','w');

k = convhull(v,w);
x1 = v(k);
y1 = w(k);
xx = floor(x1);
indexmin = find(xx==0);
[~,indexmax] = max(x1);
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),':','LineWidth',2.5); hold on;grid on

%% RS plot
t = [Rate_order1(:,1);Rate_order2(:,1)];
z = [Rate_order1(:,2);Rate_order2(:,2)];
save('t.mat','t');save('z.mat','z');

k = convhull(t,z);
x1 = t(k);
y1 = z(k);
xx = floor(x1);
indexmin = find(xx==0);
[~,indexmax] = max(x1);
plot(x1(indexmax(1):indexmin(1)),y1(indexmax(1):indexmin(1)),'*-','LineWidth',2.5);
legend('MU-LP','SC-SIC','RS');
xlabel('{\it{R_{total,1}}} (bits/s/Hz)');
ylabel('{\it{R_{total,2}}} (bits/s/Hz)');

