% RSMA CoMP JT SNR
% Using CVX for optimization (download from: cvxr.com)
function [WMMSE,p_1,p_2,p_3] = SDMA_CVX_Optimization(H,SNR,weight,g,u_opt,Rth)
u1 = weight(1);
u2 = weight(2);
u3 = weight(3);

[Nr,N_bs,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

g_1 = g(1);
g_2 = g(2);
g_3 = g(3);

u_opt_1 = u_opt(1);
u_opt_2 = u_opt(2);
u_opt_3 = u_opt(3);


%% CVX
cvx_begin quiet
variable p_1(N_bs,Nr) complex
variable p_2(N_bs,Nr) complex
variable p_3(N_bs,Nr) complex


%% objective function (i.e. WSR)
% receive power
T_1 = square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
T_2 = square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
T_3 = square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;

% eq29 MSE (function of power)
MSE_1 = abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
MSE_2 = abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;
MSE_3 = abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;

% eq33 augmented WMMSE (function of MSE)
WMMSE_1 = u_opt_1*MSE_1-log2(u_opt_1);
WMMSE_2 = u_opt_2*MSE_2-log2(u_opt_2);
WMMSE_3 = u_opt_3*MSE_3-log2(u_opt_3);

% WMMSE 
WMMSE = u1*WMMSE_1+u2*WMMSE_2+u3*WMMSE_3;

%objective function
minimize(WMMSE)

%% constraints
subject to
square_abs(p_1(1))+square_abs(p_2(1))+square_abs(p_3(1))-SNR/N_bs <= 0;
square_abs(p_1(2))+square_abs(p_2(2))+square_abs(p_3(2))-SNR/N_bs <= 0;
square_abs(p_1(3))+square_abs(p_2(3))+square_abs(p_3(3))-SNR/N_bs <= 0;
WMMSE_1-1+Rth<= 0;
WMMSE_2-1+Rth<= 0;
WMMSE_3-1+Rth<= 0;

cvx_end

end