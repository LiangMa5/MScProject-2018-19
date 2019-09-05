% RSMA SNR
% Using CVX for optimization (download from: cvxr.com)
function [WMMSE,p_1,p_2,p_3,p_12,p_13,p_23,p_123,X_1_123,X_2_123,X_3_123,X_1_12,X_2_12,X_1_13,X_3_13,X_2_23,X_3_23]=RS_CVX_Optimization(H,SNR,weight,g,u_opt,Rth)
u1 = weight(1);
u2 = weight(2);
u3 = weight(3);

[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

g_1_123 = g(1);
g_1_12 = g(2);
g_1_13 = g(3);
g_1 = g(4);

g_2_123 = g(5);
g_2_12 = g(6);
g_2_23 = g(7);
g_2 = g(8);

g_3_123 = g(9);
g_3_13 = g(10);
g_3_23 = g(11);
g_3 = g(12);

u_opt_1_123 = u_opt(1);
u_opt_1_12 = u_opt(2);
u_opt_1_13 = u_opt(3);
u_opt_1 = u_opt(4);

u_opt_2_123 = u_opt(5);
u_opt_2_12 = u_opt(6);
u_opt_2_23 = u_opt(7);
u_opt_2 = u_opt(8);

u_opt_3_123 = u_opt(9);
u_opt_3_13 = u_opt(10);
u_opt_3_23 = u_opt(11);
u_opt_3 = u_opt(12);


%% CVX
cvx_begin quiet

variable p_1(Nt,Nr) complex
variable p_2(Nt,Nr) complex
variable p_3(Nt,Nr) complex
variable p_12(Nt,Nr) complex
variable p_13(Nt,Nr) complex
variable p_23(Nt,Nr) complex
variable p_123(Nt,Nr) complex
variable X_1_123
variable X_2_123
variable X_3_123
variable X_1_12
variable X_2_12
variable X_1_13
variable X_3_13
variable X_2_23
variable X_3_23
expression constraints(1,22);

%% objective function (i.e. WMMSE)
% receive power
T_1 = square_abs(h1*p_23)+square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
T_2 = square_abs(h2*p_13)+square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
T_3 = square_abs(h3*p_12)+square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;

T_1_13 = square_abs(h1*p_13)+T_1;
T_1_12 = square_abs(h1*p_12)+T_1_13;  %s_12 is decoded before s_13
T_1_123 = square_abs(h1*p_123)+T_1_12;

T_2_23 = square_abs(h2*p_23)+T_2;
T_2_12 = square_abs(h2*p_12)+T_2_23;
T_2_123 = square_abs(h2*p_123)+T_2_12;

T_3_23 = square_abs(h3*p_23)+T_3;
T_3_13 = square_abs(h3*p_13)+T_3_23;
T_3_123 = square_abs(h3*p_123)+T_3_13;

% eq29 MSE (function of power)
MSE_1_123 = abs(g_1_123)^2*T_1_123-2*real(g_1_123*h1*p_123)+1;
MSE_1_12 = abs(g_1_12)^2*T_1_12-2*real(g_1_12*h1*p_12)+1;
MSE_1_13 = abs(g_1_13)^2*T_1_13-2*real(g_1_13*h1*p_13)+1;
MSE_1 = abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;

MSE_2_123 = abs(g_2_123)^2*T_2_123-2*real(g_2_123*h2*p_123)+1;
MSE_2_12 = abs(g_2_12)^2*T_2_12-2*real(g_2_12*h2*p_12)+1;
MSE_2_23 = abs(g_2_23)^2*T_2_23-2*real(g_2_23*h2*p_23)+1;
MSE_2 = abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;

MSE_3_123 = abs(g_3_123)^2*T_3_123-2*real(g_3_123*h3*p_123)+1;
MSE_3_13 = abs(g_3_13)^2*T_3_13-2*real(g_3_13*h3*p_13)+1;
MSE_3_23 = abs(g_3_23)^2*T_3_23-2*real(g_3_23*h3*p_23)+1;
MSE_3 = abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;

% eq33 augmented WMMSE (function of MSE)
WMMSE_1_123 = u_opt_1_123*MSE_1_123-log2(u_opt_1_123);
WMMSE_1_12 = u_opt_1_12*MSE_1_12-log2(u_opt_1_12);
WMMSE_1_13 = u_opt_1_13*MSE_1_13-log2(u_opt_1_13);
WMMSE_1 = u_opt_1*MSE_1-log2(u_opt_1);

WMMSE_2_123 = u_opt_2_123*MSE_2_123-log2(u_opt_2_123);
WMMSE_2_12 = u_opt_2_12*MSE_2_12-log2(u_opt_2_12);
WMMSE_2_23 = u_opt_2_23*MSE_2_23-log2(u_opt_2_23);
WMMSE_2 = u_opt_2*MSE_2-log2(u_opt_2);

WMMSE_3_123 = u_opt_3_123*MSE_3_123-log2(u_opt_3_123);
WMMSE_3_13 = u_opt_3_13*MSE_3_13-log2(u_opt_3_13);
WMMSE_3_23 = u_opt_3_23*MSE_3_23-log2(u_opt_3_23);
WMMSE_3 = u_opt_3*MSE_3-log2(u_opt_3);

% WMMSE (eq36 and below)
WMMSE_u1 = X_1_123+X_1_12+X_1_13+WMMSE_1;
WMMSE_u2 = X_2_123+X_2_12+X_2_23+WMMSE_2;
WMMSE_u3 = X_3_123+X_3_13+X_3_23+WMMSE_3;
WMMSE = u1*(WMMSE_u1)+u2*(WMMSE_u2)+u3*(WMMSE_u3);

%objective function
minimize(WMMSE)

%% constraints
constraints(1) = WMMSE_1_123-1-X_1_123-X_2_123-X_3_123;
constraints(2) = WMMSE_2_123-1-X_1_123-X_2_123-X_3_123;
constraints(3) = WMMSE_3_123-1-X_1_123-X_2_123-X_3_123;

constraints(4) = WMMSE_1_12-1-X_1_12-X_2_12;
constraints(5) = WMMSE_2_12-1-X_1_12-X_2_12;
constraints(6) = WMMSE_1_13-1-X_1_13-X_3_13;
constraints(7) = WMMSE_3_13-1-X_1_13-X_3_13;
constraints(8) = WMMSE_2_23-1-X_2_23-X_3_23;
constraints(9) = WMMSE_3_23-1-X_2_23-X_3_23;

constraints(10) = X_1_123;
constraints(11) = X_2_123;
constraints(12) = X_3_123;
constraints(13) = X_1_12;
constraints(14) = X_2_12;
constraints(15) = X_1_13;
constraints(16) = X_3_13;
constraints(17) = X_2_23;
constraints(18) = X_3_23;
constraints(19) = trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_3'*p_3)+trace(p_12'*p_12)+trace(p_13'*p_13)+trace(p_23'*p_23)+trace(p_123'*p_123)-SNR;
constraints(20) = WMMSE_u1-1+Rth;
constraints(21) = WMMSE_u2-1+Rth;
constraints(22) = WMMSE_u3-1+Rth;


subject to
constraints <= zeros(1,22)

cvx_end

end