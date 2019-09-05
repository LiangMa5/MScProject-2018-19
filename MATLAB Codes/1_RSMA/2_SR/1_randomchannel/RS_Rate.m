% RSMA SNR
function WSR = RS_Rate(H,SNRdB,weight,tolerance,Rth)

u1 = weight(1);
u2 = weight(2);
u3 = weight(3);

h1 = H(:,:,1);
h2 = H(:,:,2);
h3 = H(:,:,3);
    

% Decoding order of messages at layer 2: 12-->13-->23
WSR_order(1) = RS_Rate_oneorder(H,SNRdB,weight,tolerance,Rth);

% Decoding order of messages at layer 2: 12-->23-->13 
weight_2(1) = u2;
weight_2(2) = u1;
weight_2(3) = u3;
H_2(:,:,1) = h2;
H_2(:,:,2) = h1;
H_2(:,:,3) = h3;
WSR_order(2) = RS_Rate_oneorder(H_2,SNRdB,weight_2,tolerance,Rth);

% Decoding order of messages at layer 2: 13-->12-->23 
weight_3(1) = u1;
weight_3(2) = u3;
weight_3(3) = u2;
H_3(:,:,1) = h1;
H_3(:,:,2) = h3;
H_3(:,:,3) = h2;
WSR_order(3) = RS_Rate_oneorder(H_3,SNRdB,weight_3,tolerance,Rth);

% Decoding order of messages at layer 2: 13-->23-->12
weight_4(1) = u3;
weight_4(2) = u1;
weight_4(3) = u2;
H_4(:,:,1) = h3;
H_4(:,:,2) = h1;
H_4(:,:,3) = h2;
WSR_order(4) = RS_Rate_oneorder(H_4,SNRdB,weight_4,tolerance,Rth);

% Decoding order of messages at layer 2: 23-->12-->13 
weight_5(1) = u2;
weight_5(2) = u3;
weight_5(3) = u1;
H_5(:,:,1) = h2;
H_5(:,:,2) = h3;
H_5(:,:,3) = h1;
WSR_order(5) = RS_Rate_oneorder(H_5,SNRdB,weight_5,tolerance,Rth);

% Decoding order of messages at layer 2: 23-->13-->12
weight_6(1) = u3;
weight_6(2) = u2;
weight_6(3) = u1;
H_6(:,:,1) = h3;
H_6(:,:,2) = h2;
H_6(:,:,3) = h1;
WSR_order(6) = RS_Rate_oneorder(H_6,SNRdB,weight_6,tolerance,Rth);

WSR = max(WSR_order);

end





