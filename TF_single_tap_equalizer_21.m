function [est_bits,x_data] = TF_single_tap_equalizer_21(N,M,M_mod,noise_var,data_grid,Y,H_tf)
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
M_bits=log2(M_mod);
N_bits_perfram = N_syms_perfram*M_bits;
 % Y=fft(Y);
%% initial time-frequency low complexity estimate
Y_tf=fft((Y*Fn')).'; % ISFFT
X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+noise_var); % single tap equalizer
X_est = ifft(X_tf.')*Fn; % SFFT

%% detector output likelihood calculations for turbo decode
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
% scatterplot(x_data)
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
