% RSMA Rate region
% Using CVX for optimization (download from: cvxr.com)

function [WMMSE,p_1,p_2,p_c]=RS_CVX_Optimization(H,SNR,weight,g_c1,g_c2,g_p1,g_p2,u_opt_c1,u_opt_c2,u_opt_p1,u_opt_p2)
u1 = weight(1);
u2 = weight(2);

[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);                                                                                                                                                                                                                                                                                                                
 
%% CVX
cvx_begin quiet
variable p_1(Nt,Nr) complex
variable p_2(Nt,Nr) complex
variable p_c(Nt,Nr) complex

%% objective function (i.e. WMMSE)
% ---  WMMSE common message  ---
% receiver power
T_c1 = square_abs(h1*p_c)+square_abs(h1*p_1)+square_abs(h1*p_2)+1; %received power at user1
T_c2 = square_abs(h2*p_c)+square_abs(h2*p_1)+square_abs(h2*p_2)+1; %received power at user2
% eq29 MSE (function of power)
MSE_c1 = square_abs(g_c1)*T_c1-2*real(g_c1*h1*p_c)+1; %MSE of commom message at user1
MSE_c2 = square_abs(g_c2)*T_c2-2*real(g_c2*h2*p_c)+1; %MSE of commom message at user2
% eq33 WMMSE (function of MSE)
c_1 = u_opt_c1*MSE_c1-log2(u_opt_c1); 
c_2 = u_opt_c2*MSE_c2-log2(u_opt_c2); 
c = max(c_1,c_2);

%we assume common rate is for user 1, therefore, c_2=0
WMMSE_c = u1*c; %commom

% ---  WMMSE Private message  ---
% receive power
T_p1 = square_abs(h1*p_1)+square_abs(h1*p_2)+1;
T_p2 = square_abs(h2*p_1)+square_abs(h2*p_2)+1;
% eq29 MSE (function of power)
MSE_p1 = abs(g_p1)^2*T_p1-2*real(g_p1*h1*p_1)+1;
MSE_p2 = abs(g_p2)^2*T_p2-2*real(g_p2*h2*p_2)+1;
% eq33 WMMSE (function of MSE)
WMMSE_p = u1*(u_opt_p1*MSE_p1-log2(u_opt_p1))+u2*(u_opt_p2*MSE_p2-log2(u_opt_p2)); %private WMMSE

% ---  WMMSE  ---
WMMSE = WMMSE_c+WMMSE_p;
minimize(WMMSE)

%% constraints
subject to
trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_c'*p_c)-SNR <= 0;

cvx_end

end