% RSMA partial CSIT
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% subscript c - common, p - private
% Basic idea: 1. update p using g and u, 2.transform WMMSE problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder

function Rate = SDMA_Rate(H_est,Hall,SNRdB,weight,tolerance,M,c,d,Nt,N_user)
SNR = 10^(SNRdB/10);
%N_user * Nt
%% Step 1: Initialise precoder
P_private = SNR;

p_p = sqrt(P_private/N_user).*H_est'./(sum(abs(H_est).^2,2).^0.5)';

WMMSE_past = 0;  count = 0;
while true
    for i_m = 1:M % M=100, average over M realizations
        H = Hall(:,:,i_m);
        %% step 2: MMSE combinor
        % received power
        T_p = sum(abs(H*p_p).^2,2)+1;
        
        dd = abs(diag(H*p_p)).^2;
        for i_u0 = 1: N_user
            I(i_u0,:) = T_p(i_u0,:)-dd(i_u0);
        end
        
        % MMSE combinor, eq30
        g_p = diag(p_p'*H')./T_p; % gp1,gp2.. numbers
        
        %% Step 3: optimum MMSE weight
        %eq29, MMSE
        MMSE_p = T_p.\I;
        
        %eq34, optimal MMSE weight
        u_opt_p(:,:,i_m) = 1./(MMSE_p);
        
        % calculate others for cvx
        e_p(:,:,i_m) = u_opt_p(:,:,i_m).*abs(g_p).^2;%1,2...numbers
        
        for i_u = 1:N_user
            f_p(:,:,i_u,i_m) = e_p(i_u,:,i_m)*(H(i_u,:)'*H(i_u,:));
            v_p(:,:,i_u,i_m) = u_opt_p(i_u,:,i_m)*H(i_u,:)'*g_p(i_u,:)';
        end
        w_p(:,:,i_m) = log2(u_opt_p(:,:,i_m));
        
    end
    
    % sample average
    u_opt_p_aveg = mean(u_opt_p,3);
    
    e_p_aveg = mean(e_p,3);
    
    f_p_aveg = mean(f_p,4);
    
    v_p_aveg = mean(v_p,4);
    
    w_p_aveg = mean(w_p,3);
    
    %% Step 4: Update precoder, get WMMSE
    [WMMSE,p_p] = SDMA_CVX_Optimization(SNR,weight,Nt,N_user,u_opt_p_aveg,e_p_aveg,f_p_aveg,v_p_aveg,w_p_aveg);
    
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
% Rate
Rate = [w_p_aveg(1),w_p_aveg(2)];

end






