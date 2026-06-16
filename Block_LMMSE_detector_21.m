function [est_bits,x_data] = Block_LMMSE_detector_2(N,M,M_mod,noise_var,data_grid,r,H_LS1,taps)
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
%number of bits per QAM symbol
M_bits=log2(M_mod);
%number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;
%received time domain blocks 
sn_block_est=zeros(M,N);
Gn=zeros(M,M);
%% block-wise LMMSE detection
for n=1:N    
    for m=1:M
        for l=1:taps
            if(m>=l)
                Gn(m,m-l+1)=H_LS1(l,m+(n-1)*M);
            end
        end
    end
    rn=r((n-1)*M+1:n*M);    
    Rn=Gn'*Gn;
    sn_block_est(:,n)=(Rn+noise_var.*eye(M))^(-1)*(Gn'*rn);
end
X_tilda_est=sn_block_est;
%% detector output
X_est=X_tilda_est*Fn;
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
% scatterplot(x_data);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);
end
