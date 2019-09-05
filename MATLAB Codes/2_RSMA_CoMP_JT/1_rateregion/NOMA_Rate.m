% RSMA JT Rate region

function [Rate_order1,Rate_order2] = NOMA_Rate(H,SNRdB,weight,tolerance)

% Order 1: The common rate is allocated to user-1 only
Rate_order1 = NOMA_Rate_oneorder(H,SNRdB,weight,tolerance);

% Order 2: The common rate is allocated to user-2 only
H2(:,:,1) = H(:,:,2);
H2(:,:,2) = H(:,:,1);
weight2(1) = weight(2);
weight2(2) = weight(1);

Rate2 = NOMA_Rate_oneorder(H2,SNRdB,weight2,tolerance);
Rate_order2(1) = Rate2(2);
Rate_order2(2) = Rate2(1);
end
