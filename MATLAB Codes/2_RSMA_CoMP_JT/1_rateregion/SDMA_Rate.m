% RSMA CoMP JT
% simulate SDMA capacity using WMMSE Algorithm
% Implemented algorithm in the programme is from the paper:
% Weighted Sum-Rate Maximization using Weighted MMSE for MIMO-BC Beamforming Design

% notation:
% B: precoder /  A: MMSE receiver filter / W: weight matrix
% Basic idea: 1.update B using A and W, 2.transform WSR problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved based on eq21 --> get the optimum precoder
% (WSR gradient = WMMMSE gradient under eq21)
function Rk = SDMA_Capacity(H,SNRdB,weight,tolerance)

SNR = 10^(SNRdB/10);
[Nr,N_bs,N_user] = size(H);

% power for each user, WSR optimum is reached at maximum transmit power
Power_user = SNR/N_user; %Power_user(i1)
Power_bs = SNR/N_bs;

% initialise a precoder B using MRT
for i1 = 1:N_user
    h0 = H(:,:,i1);
    B0(:,:,i1) = sqrt(Power_user)*h0'/norm(h0); %%(Nt*Nbs) * Nr *Nuser  
end
B_cat = cat(2,B0(:,:,1),B0(:,:,2)); %B(1,;) --> BS1

% power constraint
Power = B_cat*B_cat';
P_bs1 = Power(1,1); % power of bs1
P_bs2 = Power(2,2);

B_catnew = [sqrt(Power_bs/P_bs1)*B_cat(1,:);sqrt(Power_bs/P_bs2)*B_cat(2,:)];

for i1 = 1:N_user
    B(:,:,i1) = B_catnew(:,i1); %(Nt*Nbs) * Nr *Nuser
end

WSR_past=100000; count=0;
while(true)
    %% Step 1: compute A(n) given B(n-1) using eq7
    for i2 = 1 : N_user
        % eq6, calculate Rvv
        Rvv(:,:,i2) = eye(Nr); % noise
        for i3 = 1:N_user
            if i3 ~= i2
                Rvv(:,:,i2) = Rvv(:,:,i2)+H(:,:,i2)*B(:,:,i3)*B(:,:,i3)'*H(:,:,i2)'; % plus interference
            end
        end
        %Rvv
        % using eq6 to get eq7
        A_mmse{i2}=B(:,:,i2)'*H(:,:,i2)'*inv(H(:,:,i2)*B(:,:,i2)*B(:,:,i2)'*H(:,:,i2)'+Rvv(:,:,i2));
    end
    
    %% Step 2: compute W(n) given B(n-1) using eq21,8
    for i4 = 1:N_user
        Ek(:,:,i4)=inv(eye(Nr)+B(:,:,i4)'*H(:,:,i4)'*inv(Rvv(:,:,i4))*H(:,:,i4)*B(:,:,i4));
        W{i4}=weight(i4)/(Ek(:,:,i4));
    end
    
    %% step 3: compute B(n) given A(n), W(n) using eq22,23
    W_blkdiag = blkdiag(W{:});
    A_mmse_blkdiag = blkdiag(A_mmse{:});
    H_all = cat(1,H(:,:,1),H(:,:,2));
    
    B_bar = inv(H_all'*A_mmse_blkdiag'*W_blkdiag*A_mmse_blkdiag*H_all+trace(W_blkdiag*A_mmse_blkdiag*A_mmse_blkdiag')*eye(N_bs)/SNR)*H_all'*A_mmse_blkdiag'*W_blkdiag;
    b = sqrt(SNR/trace(B_bar*B_bar'));
    B_mmse0 = b*B_bar; % power constraint of total power
    
    
    % power constraint of each BS
    Power = B_mmse0*B_mmse0';
    P_bs1 = Power(1,1);%norm(B_bar(1,:),'fro').^2;%Power(1,1); % power of bs1
    P_bs2 = Power(2,2);%norm(B_bar(2,:),'fro').^2;%Power(2,2);
    
    
    B_mmse = [sqrt(Power_bs/P_bs1)*B_mmse0(1,:);sqrt(Power_bs/P_bs2)*B_mmse0(2,:)];
    
    for i5 = 1:N_user
        B(:,:,i5) = B_mmse(:,(i5-1)+Nr:i5*Nr);
    end
    
    %% get capacity based on updated B using eq6,8,9
    for i6 = 1:N_user
        Rvv(:,:,i6) = eye(Nr);
        for i7 = 1:N_user
            if i7 ~= i6
                Rvv(:,:,i6)=Rvv(:,:,i6)+H(:,:,i6)*B(:,:,i7)*B(:,:,i7)'*H(:,:,i6)';
            end
        end
        Ek(:,:,i6)=inv(eye(Nr)+B(:,:,i6)'*H(:,:,i6)'/(Rvv(:,:,i6))*H(:,:,i6)*B(:,:,i6));
        Rk(i6)=real(log2(det(inv(Ek(:,:,i6)))));
    end
    
    %% compute WSR
    WSR = 0;
    for i7 = 1:N_user
        WSR = WSR+weight(i7)*Rk(i7);
    end
    
    % convergence
    if abs(WSR-WSR_past) <= tolerance
        break;
    else
        WSR_past = WSR;
        count = count+1;
    end
    
    if count >= 2000
        break;
    end
    
    
end

end