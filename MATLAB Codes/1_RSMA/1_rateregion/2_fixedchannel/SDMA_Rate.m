% RSMA Rate region
% simulate SDMA capacity using WMMSE Algorithm
% Implemented algorithm in the programme is from the paper:
% Weighted Sum-Rate Maximization using Weighted MMSE for MIMO-BC Beamforming Design

% notation:
% B: precoder /  A: MMSE receiver filter / W: weight matrix
% Basic idea: 1.update B using A and W, 2.transform WSR problem into WMMMSE problem
% the convex optimization problem (WMMSE) is solved based on eq21 --> get the optimum precoder 
% (WSR gradient = WMMMSE gradient under eq21)

function Rk = SDMA_Rate(H,SNRdB,weight,tolerance)

SNR = 10^(SNRdB/10);
% Nt is number of Tx antenna, Nr is number of Rx antenna, N_user is user number
[Nr,Nt,N_user] = size(H);

% power for each user, WSR optimum is reached at maximum transmit power
Power_user = SNR/N_user;

% initialise a precoder B using MRT
for i1 = 1:N_user
    h0 = H(:,:,i1);
    B(:,:,i1) = sqrt(Power_user)*h0'/norm(h0); %Nt * Nr *Nuser
end


WSR_past=0; count=0;
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
    
    B_bar = inv(H_all'*A_mmse_blkdiag'*W_blkdiag*A_mmse_blkdiag*H_all+trace(W_blkdiag*A_mmse_blkdiag*A_mmse_blkdiag')*eye(Nt)/SNR)*H_all'*A_mmse_blkdiag'*W_blkdiag;
    %norm(B_bar,'fro').^2-trace(B_bar*B_bar')
    b = sqrt(SNR/trace(B_bar*B_bar'));
    B_mmse = b*B_bar;   %size: Nt*(Nr*N_user)

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
    if abs(WSR-WSR_past) < tolerance
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