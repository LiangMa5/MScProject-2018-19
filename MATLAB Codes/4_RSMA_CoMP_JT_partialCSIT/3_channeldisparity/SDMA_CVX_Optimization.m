% RSMA CoMP JT partial CSIT channel disparity
% Using CVX for optimization (download from: cvxr.com)

function [WMMSE,p_p] = SDMA_CVX_Optimization(SNR,weight,N_bs,N_user,u_opt_p_aveg,e_p_aveg,f_p_aveg,v_p_aveg,w_p_aveg,Rth)

%% CVX
cvx_begin quiet
variable p_p(N_bs,N_user) complex


%% objective function (i.e. WMMSE)

% ---  WMMSE Private message  ---
for i_u = 1:N_user
    for i_sum = 1:N_user
        sum_p(i_sum) = p_p(:,i_sum)' * 0.5*(f_p_aveg(:,:,i_u)+f_p_aveg(:,:,i_u)') * p_p(:,i_sum);
    end
    WMMSE_p0(i_u,:) = sum(sum_p) + e_p_aveg(i_u) -2*real(v_p_aveg(:,:,i_u)'*p_p(:,i_u))+u_opt_p_aveg(i_u)-w_p_aveg(i_u);  
end
WMMSE_p = weight*WMMSE_p0; %private WMMSE

% ---  WMMSE  ---
WMMSE = WMMSE_p;
minimize(WMMSE)


%% constraiN_bss
subject to

for i0 = 1:N_bs
    WMMSE_p0(i0,:) <= 1-Rth;
    sum_square_abs(p_p(i0,:))-SNR/N_bs <= 0;
end

cvx_end


end