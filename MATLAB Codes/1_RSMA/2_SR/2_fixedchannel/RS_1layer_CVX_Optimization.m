% RSMA SNR
% Using CVX for optimization (download from: cvxr.com)
function [WMMSE,p_1,p_2,p_3,p_123,X_1_123,X_2_123,X_3_123]=RS_1layer_CVX_Optimization(H,SNR,weight,g,u_opt,Rth)
u1 = weight(1);
u2 = weight(2);
u3 = weight(3);

[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

g_1_123 = g(1);
g_1 = g(2);
g_2_123 = g(3);
g_2 = g(4);
g_3_123 = g(5);
g_3 = g(6);

u_opt_1_123 = u_opt(1);
u_opt_1 = u_opt(2);
u_opt_2_123 = u_opt(3);
u_opt_2 = u_opt(4);
u_opt_3_123 = u_opt(5);
u_opt_3 = u_opt(6);

%% CVX
cvx_begin quiet

variable p_1(Nt,Nr) complex
variable p_2(Nt,Nr) complex
variable p_3(Nt,Nr) complex
variable p_123(Nt,Nr) complex
variable X_1_123
variable X_2_123
variable X_3_123

expression constraints(1,10);

%% objective function (i.e. WMMSE)
% receive power
T_1 = square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
T_2 = square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
T_3 = square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;

T_1_123 = square_abs(h1*p_123)+T_1;
T_2_123 = square_abs(h2*p_123)+T_2;
T_3_123 = square_abs(h3*p_123)+T_3;

% eq29 MSE (function of power)
MSE_1_123 = abs(g_1_123)^2*T_1_123-2*real(g_1_123*h1*p_123)+1;
MSE_1 = abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;

MSE_2_123 = abs(g_2_123)^2*T_2_123-2*real(g_2_123*h2*p_123)+1;
MSE_2 = abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;

MSE_3_123 = abs(g_3_123)^2*T_3_123-2*real(g_3_123*h3*p_123)+1;
MSE_3 = abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;

% eq33 augmented WMMSE (function of MSE)
WMMSE_1_123 = u_opt_1_123*MSE_1_123-log2(u_opt_1_123);
WMMSE_1 = u_opt_1*MSE_1-log2(u_opt_1);

WMMSE_2_123 = u_opt_2_123*MSE_2_123-log2(u_opt_2_123);
WMMSE_2 = u_opt_2*MSE_2-log2(u_opt_2);

WMMSE_3_123 = u_opt_3_123*MSE_3_123-log2(u_opt_3_123);
WMMSE_3 = u_opt_3*MSE_3-log2(u_opt_3);

% WMMSE (eq36 and below)
WMMSE_u1 = X_1_123+WMMSE_1;
WMMSE_u2 = X_2_123+WMMSE_2;
WMMSE_u3 = X_3_123+WMMSE_3;
WMMSE = u1*(WMMSE_u1)+u2*(WMMSE_u2)+u3*(WMMSE_u3);

%objective function
minimize(WMMSE)

%% constraints
constraints(1) = WMMSE_1_123-1-X_1_123-X_2_123-X_3_123;
constraints(2) = WMMSE_2_123-1-X_1_123-X_2_123-X_3_123;
constraints(3) = WMMSE_3_123-1-X_1_123-X_2_123-X_3_123;

constraints(4) = X_1_123;
constraints(5) = X_2_123;
constraints(6) = X_3_123;

constraints(7) = trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_3'*p_3)+trace(p_123'*p_123)-SNR;

constraints(8) = WMMSE_u1-1+Rth;
constraints(9) = WMMSE_u2-1+Rth;
constraints(10) = WMMSE_u3-1+Rth;

subject to
constraints <= zeros(1,10)

cvx_end

end