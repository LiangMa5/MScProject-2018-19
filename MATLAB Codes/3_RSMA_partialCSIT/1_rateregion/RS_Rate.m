% RSMA CSIT Rate region

function [Rate_order1,Rate_order2] = RS_Rate(H_est,H,SNRdB,weight,tolerance,M,c,d,Nt,N_user)

% Order 1: The common rate is allocated to user-1 only
Rate_order1 = RS_1layer_Rate(H_est,H,SNRdB,weight,tolerance,M,c,d,Nt,N_user);


% Order 2: The common rate is allocated to user-2 only
weight2(1) = weight(2);
weight2(2) = weight(1);

H_est2(2,:) = H_est(1,:);
H_est2(1,:) = H_est(2,:);

for i0 = 1:M
    H2(2,:,i0) = H(1,:,i0);
    H2(1,:,i0) = H(2,:,i0);
end

Rate2 = RS_1layer_Rate(H_est2,H2,SNRdB,weight2,tolerance,M,c,d,Nt,N_user);
Rate_order2(1) = Rate2(2);
Rate_order2(2) = Rate2(1);
end
