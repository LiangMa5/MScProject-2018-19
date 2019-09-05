% RSMA SNR
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% Basic idea: 1. update p using g and u, 2.transform WSR problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder
function Rate = RS_Rate_oneorder(H,SNRdB,weight,tolerance,Rth)
SNR = 10.^(SNRdB/10);
[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

%% Step 1: Initialise precoder
% split power for 3 users
a1=0.1; a2=0.3; a3=0.6;

% precoder for 3-order stream
H_123 = [h1;h2;h3];
[~,~,V_123] = svd(H_123);  %SVD
p_123 = V_123(:,1)*sqrt(SNR*a3);

%  precoder for 2-order stream
H_12 = [h1;h2];
[~,~,V_12] = svd(H_12); %SVD
p_12 = V_12(:,1)*sqrt(SNR*a2/N_user);

H_13 = [h1;h3];
[~,~,V_13] = svd(H_13);
p_13 = V_13(:,1)*sqrt(SNR*a2/N_user);

H_23 = [h2;h3];
[~,~,V_23] = svd(H_23);
p_23 = V_23(:,1)*sqrt(SNR*a2/N_user);

% precoder for 3-order stream
p_1 = sqrt(SNR*a1/N_user)*h1'/norm(h1); %MRT
p_2 = sqrt(SNR*a1/N_user)*h2'/norm(h2);
p_3 = sqrt(SNR*a1/N_user)*h3'/norm(h3);  % size of p_i is Nt*1


WMMSE_past = 0; count = 0;
% decoding order: 12-->13-->23
while true
    %% step 2: MMSE combinor
    % received power
    T_1 = abs(h1*p_23)^2+abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
    T_2  =abs(h2*p_13)^2+abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    T_3 = abs(h3*p_12)^2+abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;
    
    T_1_13 = abs(h1*p_13)^2+T_1;
    T_1_12 = abs(h1*p_12)^2+T_1_13;  %s_12 is decoded before s_13
    T_1_123 = abs(h1*p_123)^2+T_1_12;
    
    T_2_23 = abs(h2*p_23)^2+T_2;
    T_2_12 = abs(h2*p_12)^2+T_2_23;
    T_2_123 = abs(h2*p_123)^2+T_2_12;
    
    T_3_23 = abs(h3*p_23)^2+T_3;
    T_3_13 = abs(h3*p_13)^2+T_3_23;
    T_3_123 = abs(h3*p_123)^2+T_3_13;
    
    
    % MMSE combinor, eq30
    g_1_123 = p_123'*h1'/T_1_123;
    g_1_12 = p_12'*h1'/T_1_12;
    g_1_13 = p_13'*h1'/T_1_13;
    g_1 = p_1'*h1'/T_1;
    
    g_2_123 = p_123'*h2'/T_2_123;
    g_2_12 = p_12'*h2'/T_2_12;
    g_2_23 = p_23'*h2'/T_2_23;
    g_2 = p_2'*h2'/T_2;
    
    g_3_123 = p_123'*h3'/T_3_123;
    g_3_13 = p_13'*h3'/T_3_13;
    g_3_23 = p_23'*h3'/T_3_23;
    g_3 = p_3'*h3'/T_3;
    
    g(1) = g_1_123;
    g(2) = g_1_12;
    g(3) = g_1_13;
    g(4) = g_1;
    
    g(5) = g_2_123;
    g(6) = g_2_12;
    g(7) = g_2_23;
    g(8) = g_2;
    
    g(9) = g_3_123;
    g(10) = g_3_13;
    g(11) = g_3_23;
    g(12) = g_3;
    
    %% Step 3: optimum MMSE weight
    %eq31, MMSE
    MMSE_1_123 = T_1_123\(T_1_123-abs(h1*p_123)^2);
    MMSE_1_12 = T_1_12\(T_1_12-abs(h1*p_12)^2);
    MMSE_1_13 = T_1_13\(T_1_13-abs(h1*p_13)^2);
    MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);
    
    MMSE_2_123 = T_2_123\(T_2_123-abs(h2*p_123)^2);
    MMSE_2_12 = T_2_12\(T_2_12-abs(h2*p_12)^2);
    MMSE_2_23 = T_2_23\(T_2_23-abs(h2*p_23)^2);
    MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);
    
    MMSE_3_123 = T_3_123\(T_3_123-abs(h3*p_123)^2);
    MMSE_3_13 = T_3_13\(T_3_13-abs(h3*p_13)^2);
    MMSE_3_23 = T_3_23\(T_3_23-abs(h3*p_23)^2);
    MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);
    
    %eq34, optimal MMSE weight
    u_opt_1_123 = inv(MMSE_1_123);
    u_opt_1_12 = inv(MMSE_1_12);
    u_opt_1_13 = inv(MMSE_1_13);
    u_opt_1 = inv(MMSE_1);
    
    u_opt_2_123 = inv(MMSE_2_123);
    u_opt_2_12 = inv(MMSE_2_12);
    u_opt_2_23 = inv(MMSE_2_23);
    u_opt_2 = inv(MMSE_2);
    
    u_opt_3_123 = inv(MMSE_3_123);
    u_opt_3_13 = inv(MMSE_3_13);
    u_opt_3_23 = inv(MMSE_3_23);
    u_opt_3 = inv(MMSE_3);
    
    u_opt(1) = u_opt_1_123;
    u_opt(2) = u_opt_1_12;
    u_opt(3) = u_opt_1_13;
    u_opt(4) = u_opt_1;
    
    u_opt(5) = u_opt_2_123;
    u_opt(6) = u_opt_2_12;
    u_opt(7) = u_opt_2_23;
    u_opt(8) = u_opt_2;
    
    u_opt(9) = u_opt_3_123;
    u_opt(10) = u_opt_3_13;
    u_opt(11) = u_opt_3_23;
    u_opt(12) = u_opt_3;
    
    
    %% Step 4: Update precoder, get WSR
    [WMMSE,p_1,p_2,p_3,p_12,p_13,p_23,p_123,X_1_123,X_2_123,X_3_123,X_1_12,X_2_12,X_1_13,X_3_13,X_2_23,X_3_23]=RS_CVX_Optimization(H,SNR,weight,g,u_opt,Rth);
    if abs(WMMSE-WMMSE_past) <= tolerance
        break;
    else
        WMMSE_past = WMMSE;
        count = count+1;
    end
    
    if count >= 2000
        break;
    end
    
end


%% Calculate the rate of each UE
%Private rate 
T_1 = abs(h1*p_23)^2+abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
T_2 = abs(h2*p_13)^2+abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
T_3 = abs(h3*p_12)^2+abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;

MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);
MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);
MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);

R_1 = -real(log2(MMSE_1));
R_2 = -real(log2(MMSE_2));
R_3 = -real(log2(MMSE_3));

%Common rate 
R_1_123=-X_1_123;
R_2_123=-X_2_123;
R_3_123=-X_3_123;
R_1_12 =-X_1_12;
R_2_12 =-X_2_12;
R_1_13 =-X_1_13;
R_3_13 =-X_3_13;
R_2_23 =-X_2_23;
R_3_23 =-X_3_23;

% Rate of each user
c_1 = R_1_123 + R_1_12 + R_1_13 + R_1;
c_2 = R_2_123 + R_2_12 + R_2_23 + R_2;
c_3 = R_3_123 + R_3_13 + R_3_23 + R_3;

% WSR
Rate = weight(1)*c_1+weight(2)*c_2+weight(3)*c_3;

end







