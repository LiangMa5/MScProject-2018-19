% RSMA SNR
% Using CVX for optimization (download from: cvxr.com)
function [WMMSE,p_1,p_2,p_3]=NOMA_CVX_Optimization(H,SNR,weight,g,u_opt,Rth)
u1 = weight(1);
u2 = weight(2);
u3 = weight(3);

[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

g_1 = g(1);
g_2_1 = g(2);
g_2 = g(3);
g_3_1 = g(4);
g_3_2 = g(5);
g_3 = g(6);

u_opt_1 = u_opt(1);
u_opt_2_1 = u_opt(2);
u_opt_2 = u_opt(3);
u_opt_3_1 = u_opt(4);
u_opt_3_2 = u_opt(5);
u_opt_3 = u_opt(6);


%% CVX
cvx_begin quiet

variable p_1(Nt,Nr) complex
variable p_2(Nt,Nr) complex
variable p_3(Nt,Nr) complex

%% objective function (i.e. WMMSE)
% receive power
T_1 = square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;

T_2_1 = square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
T_2 = square_abs(h2*p_2)+square_abs(h2*p_3)+1;

T_3_1 = square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;
T_3_2 = square_abs(h3*p_2)+square_abs(h3*p_3)+1;
T_3 = square_abs(h3*p_3)+1;

% eq29 MSE (function of power)
MMSE_1 = abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
MMSE_2_1 = abs(g_2_1)^2*T_2_1-2*real(g_2_1*h2*p_1)+1;
MMSE_2 = abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;

MMSE_3_1 = abs(g_3_1)^2*T_3_1-2*real(g_3_1*h3*p_1)+1;
MMSE_3_2 = abs(g_3_2)^2*T_3_2-2*real(g_3_2*h3*p_2)+1;
MMSE_3 = abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;

% eq33 augmented WMMSE (function of MSE)
WMMSE_1 = u_opt_1*MMSE_1-log2(u_opt_1);
WMMSE_2_1 = u_opt_2_1*MMSE_2_1-log2(u_opt_2_1);
WMMSE_2 = u_opt_2*MMSE_2-log2(u_opt_2);

WMMSE_3_1 = u_opt_3_1*MMSE_3_1-log2(u_opt_3_1);
WMMSE_3_2 = u_opt_3_2*MMSE_3_2-log2(u_opt_3_2);
WMMSE_3 = u_opt_3*MMSE_3-log2(u_opt_3);


% WMMSE 
WMMSE_u1 = max([WMMSE_1,WMMSE_2_1,WMMSE_3_1]);
WMMSE_u2 = max([WMMSE_2,WMMSE_3_2]);
WMMSE_u3 = WMMSE_3;
WMMSE = u1*(WMMSE_u1)+u2*(WMMSE_u2)+u3*(WMMSE_u3);

%objective function
minimize(WMMSE)

%% constraints
subject to
trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_3'*p_3)-SNR <= 0

WMMSE_u1-1+Rth<= 0;
WMMSE_u2-1+Rth<= 0;
WMMSE_u3-1+Rth<= 0;

cvx_end

end