function [est_bits,ite,x_data] = MRC_delay_time_detector_21(N,M,M_data,M_mod,no,data_grid,r,H_tf,nu_ml_tilda,taps,omega,decision,init_estimate,n_ite_MRC)
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
%% Initial assignments
% %Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
%number of bits per QAM symbol
M_bits=log2(M_mod);
%number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;
%received delay-time samples 
Y_tilda=reshape(r,M,N);
 M_prime=M_data;
%% initial estimate using single tap TF equalizer
if(init_estimate==1)
    Y_tf=fft(Y_tilda).'; % delay-time to frequency-time domain                                                                      % equation (63) in [R1]                                   
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no); % single tap equalizer                                                            % equation (64) in [R1]
    X_est = ifft(X_tf.')*Fn; % SFFT                                                                                                 % equation (65) in [R1]
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
    X_tilda_est=X_est*Fn';
else
    X_est=zeros(M,N);
    X_tilda_est=X_est*Fn';
end
x_m=X_est.';
x_m_tilda=X_tilda_est.';


%% MRC detector    %% Algorithm 2 in [R1] (or Algotithm 3 in Chapter 6, [R2])
%% initial computation
d_m_tilda=zeros(N,M);
y_m_tilda=reshape(r,M,N).';
delta_y_m_tilda=y_m_tilda;
for m=1:M_prime   
    for l=1:taps
        d_m_tilda(:,m)=d_m_tilda(:,m)+abs(nu_ml_tilda(:,m+(l-1),l).^2);                                                             % equation (59) in [R1]
    end
end
for m=1:M         
    for l=1:taps
        if(m>=l)
            delta_y_m_tilda(:,m)=delta_y_m_tilda(:,m)-nu_ml_tilda(:,m,l).*x_m_tilda(:,m-(l-1));                                     % Line 3 of Algorithm 2 in [R1]
        end
    end
end
x_m_tilda_old=x_m_tilda;
c_m_tilda=x_m_tilda;

%% iterative computation
for ite=1:n_ite_MRC                                                                                                                 % Line 5 of Algorithm 2 in [R1]
    delta_g_m_tilda=zeros(N,M);
    for m=1:M_prime         
        for l=1:taps
            delta_g_m_tilda(:,m)=delta_g_m_tilda(:,m)+conj(nu_ml_tilda(:,m+(l-1),l)).*delta_y_m_tilda(:,m+(l-1));                   % Line 8 of Algorithm 2 in [R1]
        end
        c_m_tilda(:,m)=x_m_tilda_old(:,m)+delta_g_m_tilda(:,m)./d_m_tilda(:,m);                                                     % Line 9 of Algorithm 2 in [R1]
        if(decision==1)
            x_m(:,m)=qammod(qamdemod(Fn*(c_m_tilda(:,m)),M_mod,'gray'),M_mod,'gray');                                               % Line 10 of Algorithm 2 in [R1]
            x_m_tilda(:,m)=(1-omega)*c_m_tilda(:,m)+omega*Fn'*x_m(:,m);
        else
            x_m_tilda(:,m)=c_m_tilda(:,m);
        end
        for l=1:taps                                                                                                                 % Line 11 of Algorithm 2 in [R1]
            delta_y_m_tilda(:,m+(l-1))=delta_y_m_tilda(:,m+(l-1))-nu_ml_tilda(:,m+(l-1),l).*(x_m_tilda(:,m)-x_m_tilda_old(:,m));    % Line 12 of Algorithm 2 in [R1]
        end                                                                                                                         % Line 13 of Algorithm 2 in [R1] 
        x_m_tilda_old(:,m)=x_m_tilda(:,m);
    end
       
    %% convergence criteria
    errors(ite)=norm(delta_y_m_tilda);
    if(ite>1)
        if(errors(ite)>=errors(ite-1))                                                                                                 % Line 15 of Algorithm 2 in [R1]
            break;
        end
    end   
end
if(n_ite_MRC==0)
    ite=0;
end
%% detector output bits
X_est=(Fn*x_m_tilda).';
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
 % x_data=x_data.*4;
    % scatterplot(x_data)
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);                                                % Line 17 of Algorithm 2 in [R1]

end
