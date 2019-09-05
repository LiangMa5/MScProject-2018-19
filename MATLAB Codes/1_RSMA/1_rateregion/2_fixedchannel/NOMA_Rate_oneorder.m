% RSMA Rate region
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% subscript c - common, p - private
% Basic idea: 1. update p using g and u, 2.transform WMMSE problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder

function Rate = NOMA_Rate_oneorder(H,SNRdB,weight,tolerance)
SNR = 10^(SNRdB/10);
[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);


%% Step 1: Initialise precoder
P_ue = (SNR)/N_user;
p_1 = h1'/norm(h1)*sqrt(P_ue);
p_2 = h2'/norm(h2)*sqrt(P_ue);


WMMSE_past = 0;
count = 0;
while true
    %% step 2: MMSE combinor
    % received power
    T_1_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_2_1 = abs(h2*p_1)^2+abs(h2*p_2)^2+1;
    T_2_2 = abs(h2*p_2)^2+1;
    I_1_1 = abs(h1*p_2)^2+1;
    I_2_1 = abs(h2*p_2)^2+1;
    I_2_2=1;
    
    % MMSE combinor, eq30
    g_1_1 = p_1'*h1'/T_1_1;
    g_2_1 = p_1'*h2'/T_2_1;
    g_2_2 = p_2'*h2'/T_2_2;
    
    %% Step 3: optimum MMSE weight
    %eq31, MMSE
    MMSE_1_1 = T_1_1\I_1_1;
    MMSE_2_1 = T_2_1\I_2_1;
    MMSE_2_2 = T_2_2\I_2_2;
    
    %eq34, optimal MMSE weight
    u_opt_1_1 = inv(MMSE_1_1);
    u_opt_2_1 = inv(MMSE_2_1);
    u_opt_2_2 = inv(MMSE_2_2);
    
    
    %optimization
    [WMMSE,p_1,p_2] = NOMA_CVX_Optimization(weight,H,SNR,g_1_1,g_2_1,g_2_2,u_opt_1_1,u_opt_2_1,u_opt_2_2);
    
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
% receive power
T_2_2 = abs(h2*p_2)^2+1;
I_2_2 = 1;
T_1_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+1;
T_2_1 = abs(h2*p_1)^2+abs(h2*p_2)^2+1;
I_1_1 = abs(h1*p_2)^2+1;
I_2_1 = abs(h2*p_2)^2+1;


R_1_1 = real(log2(T_1_1/I_1_1));
R_2_1 = real(log2(T_2_1/I_2_1));

%Private Rate
R_2 = real(log2(T_2_2/I_2_2));

Rate(1,:) = min(R_1_1,R_2_1);
Rate(2,:) = R_2;

end







