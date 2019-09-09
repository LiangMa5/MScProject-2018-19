% RSMA CoMP JT partial CSIT
% notation: T: power / p: precoder / g: combinor / u: MMSE weight
% subscript c - common, p - private
% Basic idea: 1. update p using g and u, 2.transform WMMSE problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved using CVX --> get the optimum precoder

function Rate = RS_1layer_Rate(H_est,Hall,SNRdB,weight,tolerance,M,c,d,N_bs,N_user,Rth)
SNR = 10^(SNRdB/10);
%N_user * Nt
%% Step 1: Initialise precoder
P_common = SNR-SNR^(d);
P_private = SNR^(d);

[~,~,V] = svd(H_est); %SVD, V is the precoder (see papar section V.A)
p_c = sqrt(P_common)*V(:,1);
p_p = sqrt(P_private/N_user).*H_est'./(sum(abs(H_est).^2,2).^0.5)';

% %two user
% h1 = H_est(1,:);
% h2 = H_est(2,:);
% if norm(h1) >= norm(h2)
%     p_p = sqrt(P_private.*[0.9,0.1]).*H_est'./(sum(abs(H_est).^2,2).^0.5)';
% end
% 
% if norm(h1) < norm(h2)
%     p_p = sqrt(P_private.*[0.1,0.9]).*H_est'./(sum(abs(H_est).^2,2).^0.5)';
% end


P0 = [p_c,p_p];
Power_bs = SNR/N_bs;
Power = P0*P0';

for i0 = 1:N_bs
    P_bs(i0,:) = Power (i0,i0);
end
P = sqrt(Power_bs)*P0./sqrt(P_bs);

WMMSE_past = 0;  count = 0;
while true
    for i_m = 1:M % M=100, average over M realizations
        H = Hall(:,:,i_m);
        %% step 2: MMSE combinor
        % received power
        T_c = sum(abs(H*P).^2,2)+1; %1;2;....
        T_p = sum(abs(H*p_p).^2,2)+1;
          
        dd = abs(diag(H*p_p)).^2;
        for i_u0 = 1: N_user       
            I(i_u0,:) = T_p(i_u0,:)-dd(i_u0);
        end
        
   
        % MMSE combinor, eq30
        g_c = ((p_c'*H').')./T_c;
        g_p = diag(p_p'*H')./T_p; % gp1,gp2.. numbers
        
        
        %% Step 3: optimum MMSE weight
        %eq31, MMSE
        MMSE_c = T_c.\T_p;
        MMSE_p = T_p.\I;
        
        %eq34, optimal MMSE weight
        u_opt_c(:,:,i_m) = 1./(MMSE_c);
        u_opt_p(:,:,i_m) = 1./(MMSE_p);
        
        % calculate others for cvx
        e_c(:,:,i_m) = u_opt_c(:,:,i_m).*abs(g_c).^2;
        e_p(:,:,i_m) = u_opt_p(:,:,i_m).*abs(g_p).^2;%1,2...numbers
                
        for i_u = 1:N_user
            f_c(:,:,i_u,i_m) = e_c(i_u,:,i_m)*(H(i_u,:)'*H(i_u,:)); %matrix
            f_p(:,:,i_u,i_m) = e_p(i_u,:,i_m)*(H(i_u,:)'*H(i_u,:));
            
            v_c(:,:,i_u,i_m) = u_opt_c(i_u,:,i_m)*H(i_u,:)'*g_c(i_u,:)'; %matrix
            v_p(:,:,i_u,i_m) = u_opt_p(i_u,:,i_m)*H(i_u,:)'*g_p(i_u,:)';
        end
        
        w_c(:,:,i_m) = log2(u_opt_c(:,:,i_m));
        w_p(:,:,i_m) = log2(u_opt_p(:,:,i_m));
        
    end
    
    % sample average
    u_opt_c_aveg = mean(u_opt_c,3);
    u_opt_p_aveg = mean(u_opt_p,3);
    
    e_c_aveg = mean(e_c,3);
    e_p_aveg = mean(e_p,3);
    
    f_c_aveg = mean(f_c,4);
    f_p_aveg = mean(f_p,4);
    
    v_c_aveg = mean(v_c,4);
    v_p_aveg = mean(v_p,4);
    
    w_c_aveg = mean(w_c,3);
    w_p_aveg = mean(w_p,3);
    

    %% Step 4: Update precoder, get WMMSE
    [WMMSE,p_p,p_c] = RS_1layer_CVX_Optimization(SNR,weight,N_bs,N_user,u_opt_c_aveg,u_opt_p_aveg,e_c_aveg,e_p_aveg,f_c_aveg,f_p_aveg,v_c_aveg,v_p_aveg,w_c_aveg,w_p_aveg,Rth);
    
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
Rate = min(w_c_aveg)+ sum(w_p_aveg);


end






