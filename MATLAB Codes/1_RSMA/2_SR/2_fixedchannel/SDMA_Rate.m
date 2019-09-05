% RSMA SNR
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% Basic idea: 1. update p using g and u, 2.transform WSR problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder
function Rate = SDMA_Rate(H,SNRdB,weight,tolerance,Rth)
SNR = 10.^(SNRdB/10);
[Nr,Nt,N_user] = size(H);
h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);

%% Step 1: Initialise precoder
% precoder 
a1=0.3; a2=0.3; a3=0.4;
p_1 = sqrt(SNR*a1)*h1'/norm(h1); %MRT
p_2 = sqrt(SNR*a2)*h2'/norm(h2);
p_3 = sqrt(SNR*a3)*h3'/norm(h3);  % size of p_i is Nt*1


WMMSE_past = 0; count = 0;
while true
    %% step 2: MMSE combinor
    % received power
    T_1 = abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
    T_2 = abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    T_3 = abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;
       
    % MMSE combinor, eq30
    g_1 = p_1'*h1'/T_1;
    g_2 = p_2'*h2'/T_2;    
    g_3 = p_3'*h3'/T_3;
    
    g(1) = g_1;
    g(2) = g_2;
    g(3) = g_3;
    
    %% Step 3: optimum MMSE weight
    % eq31, MMSE  
    MMSE_1 = T_1\(T_1-abs(h1*p_1)^2);   
    MMSE_2 = T_2\(T_2-abs(h2*p_2)^2);    
    MMSE_3 = T_3\(T_3-abs(h3*p_3)^2);
    
    % eq34, optimal MMSE weight
    u_opt_1 = inv(MMSE_1);  
    u_opt_2 = inv(MMSE_2);
    u_opt_3 = inv(MMSE_3);
    
    u_opt(1) = u_opt_1;  
    u_opt(2) = u_opt_2;
    u_opt(3) = u_opt_3;
    
    
    %% Step 4: Update precoder, get WSR
    [WMMSE,p_1,p_2,p_3] = SDMA_CVX_Optimization(H,SNR,weight,g,u_opt,Rth);
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

% WSR
Rate = weight(1)*R_1+weight(2)*R_2+weight(3)*R_3;

end







