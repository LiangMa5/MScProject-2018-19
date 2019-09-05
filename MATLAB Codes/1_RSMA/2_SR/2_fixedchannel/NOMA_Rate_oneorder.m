% RSMA SNR
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% Basic idea: 1. update p using g and u, 2.transform WSR problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder
function Rate = NOMA_Rate_oneorder(H,SNRdB,weight,tolerance,Rth)
SNR = 10.^(SNRdB/10);
[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

%% Step 1: Initialise precoder
% split power for 3 users
a1=0.4; a2=0.3; a3=0.3;

% precoder for 3-order stream
% decoding order 1-->2-->3
H_123 = [h1;h2;h3];
[~,~,V_123] = svd(H_123);  %SVD
p_1 = V_123(:,1)*sqrt(SNR*a1);

H_23 = [h2;h3];
[~,~,V_23] = svd(H_23);
p_2 = V_23(:,1)*sqrt(SNR*a2);

p_3 = sqrt(SNR*a3)*h3'/norm(h3); %MRT


WMMSE_past = 0; count = 0;
% decoding order: 12-->13-->23
while true
    %% step 2: MMSE combinor
    % received power
    T_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
    
    T_2_1 = abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    T_2 = abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    
    T_3_1 = abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;
    T_3_2 = abs(h3*p_2)^2+abs(h3*p_3)^2+1;
    T_3 = abs(h3*p_3)^2+1;
    
    
    % MMSE combinor, eq30
    g_1 = p_1'*h1'/T_1;
    
    g_2_1 = p_1'*h2'/T_2_1;
    g_2 = p_2'*h2'/T_2;
    
    g_3_1 = p_1'*h3'/T_3_1;
    g_3_2 = p_2'*h3'/T_3_2;
    g_3 = p_3'*h3'/T_3;
    
    g(1) = g_1;  
    g(2) = g_2_1;
    g(3) = g_2;
    g(4) = g_3_1;
    g(5) = g_3_2;
    g(6) = g_3;
    
    %% Step 3: optimum MMSE weight
    % eq31, MMSE
    MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);
    MMSE_2_1 = T_2_1\(T_2_1-abs(h2*p_1)^2);
    MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);
    MMSE_3_1 = T_3_1\(T_3_1-abs(h3*p_1)^2);
    MMSE_3_2 = T_3_2\(T_3_2-abs(h3*p_2)^2);
    MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);
    
    % eq34, optimal MMSE weight
    u_opt_1 = inv(MMSE_1);
    u_opt_2_1 = inv(MMSE_2_1);
    u_opt_2 = inv(MMSE_2);
    u_opt_3_1 = inv(MMSE_3_1);
    u_opt_3_2 = inv(MMSE_3_2);
    u_opt_3 = inv(MMSE_3);
    
    u_opt(1) = u_opt_1;
    u_opt(2) = u_opt_2_1;
    u_opt(3) = u_opt_2;
    u_opt(4) = u_opt_3_1;
    u_opt(5) = u_opt_3_2;
    u_opt(6) = u_opt_3;
    
    
    %% Step 4: Update precoder, get WSR
    [WMMSE,p_1,p_2,p_3]=NOMA_CVX_Optimization(H,SNR,weight,g,u_opt,Rth);
    
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
% received power
T_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;

T_2_1 = abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
T_2 = abs(h2*p_2)^2+abs(h2*p_3)^2+1;

T_3_1 = abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;
T_3_2 = abs(h3*p_2)^2+abs(h3*p_3)^2+1;
T_3 = abs(h3*p_3)^2+1;



% eq31, MMSE
MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);
MMSE_2_1 = T_2_1\(T_2_1-abs(h2*p_1)^2);
MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);
MMSE_3_1 = T_3_1\(T_3_1-abs(h3*p_1)^2);
MMSE_3_2 = T_3_2\(T_3_2-abs(h3*p_2)^2);
MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);


R_1 = -real(log2(MMSE_1));
R_2_1 = -real(log2(MMSE_2_1));
R_2 = -real(log2(MMSE_2));
R_3_1 = -real(log2(MMSE_3_1));
R_3_2 = -real(log2(MMSE_3_2));
R_3 = -real(log2(MMSE_3));


c_u1 = min([R_1,R_2_1,R_3_1]);
c_u2 = min([R_2,R_3_2]);
c_u3 = R_3;
% WSR
Rate = weight(1)*c_u1+weight(2)*c_u2+weight(3)*c_u3;

end







