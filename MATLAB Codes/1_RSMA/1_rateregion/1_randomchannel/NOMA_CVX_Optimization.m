% RSMA
% Using CVX for optimization (download from: cvxr.com)

function [WMMSE,p_1,p_2] = NOMA_CVX_Optimization(weight,H,SNR,g_1_1,g_2_1,g_2_2,u_opt_1_1,u_opt_2_1,u_opt_2_2)
u1 = weight(1);
u2 = weight(2);

[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);

%% CVX
cvx_begin quiet
variable p_1(Nt,Nr) complex
variable p_2(Nt,Nr) complex

%% objective function (i.e. WMMSE)
% receiver power
T_1_1 = square_abs(h1*p_1)+square_abs(h1*p_2)+1;
T_2_1 = square_abs(h2*p_1)+square_abs(h2*p_2)+1;
% eq29 MSE (function of power)
MSE_1_1 = square_abs(g_1_1)*T_1_1-2*real(g_1_1*h1*p_1)+1;
MSE_2_1 = square_abs(g_2_1)*T_2_1-2*real(g_2_1*h2*p_1)+1;
% eq33 WMMSE (function of MSE)
c_1_1=(u_opt_1_1*MSE_1_1-log2(u_opt_1_1));
c_2_1 =(u_opt_2_1*MSE_2_1-log2(u_opt_2_1));
c = max(c_1_1,c_2_1);
WMMSE_c = u1*c;

% receive power
T_2_2 = square_abs(h2*p_2)+1;
% eq29 MSE (function of power)
MSE_2_2 = abs(g_2_2)^2*T_2_2-2*real(g_2_2*h2*p_2)+1;
WMMSE_p = u2*(u_opt_2_2*MSE_2_2-log2(u_opt_2_2));

WMMSE = WMMSE_c+WMMSE_p;
minimize(WMMSE)

%% constraints
subject to
trace(p_1'*p_1)+trace(p_2'*p_2)-SNR <= 0;

cvx_end

end