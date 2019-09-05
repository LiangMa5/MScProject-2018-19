% RSMA CoMP JT SNR
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% Basic idea: 1. update p using g and u, 2.transform WSR problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder
function Rate = RS_1layer_Rate(H,SNRdB,weight,tolerance,Rth)
SNR = 10.^(SNRdB/10);
[Nr,N_bs,N_user] = size(H);
Power_bs = SNR/N_bs;

h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

%% Step 1: Initialise precoder
% split power for 3 users
a1=0.2; a3=1-a1;

% precoder for 3-order stream
H_123 = [h1;h2;h3];
[~,~,V_123] = svd(H_123);  %SVD
p_123 = V_123(:,1)*sqrt(SNR*a3);

% precoder for 3-order stream
p_1 = sqrt(SNR*a1/N_user)*h1'/norm(h1); %MRT
p_2 = sqrt(SNR*a1/N_user)*h2'/norm(h2);
p_3 = sqrt(SNR*a1/N_user)*h3'/norm(h3);  % size of p_i is Nt*1


% power constraint
p_cat = [p_123,p_1,p_2,p_3];
Power = p_cat*p_cat';
P_bs1 = Power(1,1); % power of bs1
P_bs2 = Power(2,2);
P_bs3 = Power(3,3);
p_catnew = [sqrt(Power_bs/P_bs1)*p_cat(1,:);sqrt(Power_bs/P_bs2)*p_cat(2,:);sqrt(Power_bs/P_bs3)*p_cat(3,:)];

p_123 = p_catnew(:,1);
p_1 = p_catnew(:,1);
p_2 = p_catnew(:,2);
p_3 = p_catnew(:,3);


WMMSE_past = 0; count = 0;
while true
    %% step 2: MMSE combinor
    % received power
    T_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
    T_2 = abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    T_3 = abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;

    T_1_123 = abs(h1*p_123)^2+T_1;    
    T_2_123 = abs(h2*p_123)^2+T_2;
    T_3_123 = abs(h3*p_123)^2+T_3;
   
    % MMSE combinor, eq30
    g_1_123 = p_123'*h1'/T_1_123;
    g_1 = p_1'*h1'/T_1;
    
    g_2_123 = p_123'*h2'/T_2_123;
    g_2 = p_2'*h2'/T_2;
    
    g_3_123 = p_123'*h3'/T_3_123;
    g_3 = p_3'*h3'/T_3;
    
    g(1) = g_1_123;
    g(2) = g_1;    
    g(3) = g_2_123;
    g(4) = g_2;    
    g(5) = g_3_123;
    g(6) = g_3;
    
    %% Step 3: optimum MMSE weight
    %eq31, MMSE
    MMSE_1_123 = T_1_123\(T_1_123-abs(h1*p_123)^2);
    MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);
    
    MMSE_2_123 = T_2_123\(T_2_123-abs(h2*p_123)^2);
    MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);
    
    MMSE_3_123 = T_3_123\(T_3_123-abs(h3*p_123)^2);
    MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);
    
    %eq34, optimal MMSE weight
    u_opt_1_123 = inv(MMSE_1_123);
    u_opt_1 = inv(MMSE_1);
    
    u_opt_2_123 = inv(MMSE_2_123);
    u_opt_2 = inv(MMSE_2);
    
    u_opt_3_123 = inv(MMSE_3_123);
    u_opt_3 = inv(MMSE_3);
    
    u_opt(1) = u_opt_1_123;
    u_opt(2) = u_opt_1;   
    u_opt(3) = u_opt_2_123;
    u_opt(4) = u_opt_2;
    u_opt(5) = u_opt_3_123;
    u_opt(6) = u_opt_3;
    
    
    %% Step 4: Update precoder, get WSR
    [WMMSE,p_1,p_2,p_3,p_123,X_1_123,X_2_123,X_3_123]=RS_1layer_CVX_Optimization(H,SNR,weight,g,u_opt,Rth);
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
T_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
T_2 = abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
T_3 = abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;

MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);
MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);
MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);

R_1 = -real(log2(MMSE_1));
R_2 = -real(log2(MMSE_2));
R_3 = -real(log2(MMSE_3));

%Common rate 
R_1_123 = -X_1_123;
R_2_123 = -X_2_123;
R_3_123 = -X_3_123;

% Rate of each user
c_1 = R_1_123 + R_1;
c_2 = R_2_123 + R_2;
c_3 = R_3_123 + R_3;

% WSR
Rate = weight(1)*c_1+weight(2)*c_2+weight(3)*c_3;

end







