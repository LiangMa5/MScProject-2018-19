% RSMA CoMP JT Rate region
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% subscript c - common, p - private
% Basic idea: 1. update p using g and u, 2.transform WMMSE problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder

function Rate = RS_Rate_oneorder(H,SNRdB,weight,tolerance)
SNR = 10^(SNRdB/10);
[Nr,N_bs,N_user] = size(H);
Power_bs = SNR/N_bs;
h1 = H(:,:,1);
h2 = H(:,:,2);
H_concat = [h1;h2];

%% Step 1: Initialise precoder
a1 = 0.8;
P_common = SNR*a1;
P_private = SNR*(1-a1);
 
[~,~,V] = svd(H_concat); %SVD, V is the precoder (see papar section V.A)
p_c = sqrt(P_common)*V(:,1);

p_1 = sqrt(P_private/N_user)*h1'/norm(h1); % using MRT 
p_2 = sqrt(P_private/N_user)*h2'/norm(h2);

% power constraint
p_cat = [p_1,p_2,p_c];
Power = p_cat*p_cat';
P_bs1 = Power(1,1); % power of bs1
P_bs2 = Power(2,2);

p_catnew = [sqrt(Power_bs/P_bs1)*p_cat(1,:);sqrt(Power_bs/P_bs2)*p_cat(2,:)];

p_1 = p_catnew(:,1);
p_2 = p_catnew(:,2);
p_c = p_catnew(:,3);

WMMSE_past = 0;
count = 0;
while true
    %% step 2: MMSE combinor
    % received power
    T_c1 = abs(h1*p_c)^2+abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_c2 = abs(h2*p_c)^2+abs(h2*p_1)^2+abs(h2*p_2)^2+1;
    
    T_p1 = abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_p2 = abs(h2*p_1)^2+abs(h2*p_2)^2+1;
    I_1 = abs(h1*p_2)^2+1;
    I_2 = abs(h2*p_1)^2+1;
    
    % MMSE combinor, eq30
    g_c1 = p_c'*h1'/T_c1;
    g_c2 = p_c'*h2'/T_c2;
    
    g_p1 = p_1'*h1'/T_p1;
    g_p2 = p_2'*h2'/T_p2;
   

    %% Step 3: optimum MMSE weight
    %eq31, MMSE
    MMSE_c1 = T_c1\T_p1;
    MMSE_c2 = T_c2\T_p2;
    MMSE_p1 = T_p1\I_1;
    MMSE_p2 = T_p2\I_2;
       
    %eq34, optimal MMSE weight
    u_opt_c1 = inv(MMSE_c1);
    u_opt_c2 = inv(MMSE_c2);
    u_opt_p1 = inv(MMSE_p1);
    u_opt_p2 = inv(MMSE_p2);
    

    %% Step 4: Update precoder, get WMMSE
    [WMMSE,p_1,p_2,p_c] = RS_CVX_Optimization(weight,H,SNR,g_c1,g_c2,g_p1,g_p2,u_opt_c1,u_opt_c2,u_opt_p1,u_opt_p2);
    
    % convergence
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
T_p1 = abs(h1*p_1)^2+abs(h1*p_2)^2+1;
T_p2 = abs(h2*p_1)^2+abs(h2*p_2)^2+1;
T_c1 = abs(h1*p_c)^2+T_p1;
T_c2 = abs(h2*p_c)^2+T_p2;
I_1 = T_p1-abs(h1*p_1)^2;
I_2 = T_p2-abs(h2*p_2)^2;

% Common Rate
R_c1 = real(log2(T_c1/T_p1));
R_c2 = real(log2(T_c2/T_p2));

% Private Rate
R_p1 = real(log2(T_p1/I_1));
R_p2 = real(log2(T_p2/I_2));

% Rate
Rate(1,:) = R_p1 + min(R_c1,R_c2); %assume common stream is for user1 only
Rate(2,:) = R_p2;


end



          
          
          
 