function [est_bits,ite,x_est] = MPA_detector_21(N,M,M_mod,no,data_grid,y,H,n_ite)

%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
% data start delay index in the delay-Doppler grid
M_bits=log2(M_mod);
N_bits_perfram = N_syms_perfram*M_bits;

alphabet=qammod(0:M_mod-1,M_mod,'gray').';
x_est=zeros(M*N,1);
ind_ba=zeros(N,M*N);
ind_ab=zeros(N,M*N);
length_ba=zeros(M*N,1);
length_ab=zeros(M*N,1);
for b=1:N*M
    h_ba=H(b,:);
    ind_sparse=find(abs(h_ba)>1e-6);
    ind_sparse=sort(ind_sparse);
    length_ba(b)=length(ind_sparse);
    ind_ba(1:length_ba(b),b)=ind_sparse;
end
for a=1:N*M
    h_ab=H(:,a);
    ind_sparse=find(abs(h_ab)>1e-6);
    length_ab(a)=length(ind_sparse);
    ind_ab(1:length_ab(a),a)=ind_sparse;
end

taps=max(max(length_ab),max(length_ba));
p_map = ones(N*M,N*M,M_mod)*(1/M_mod);

delta_fra=0.5;
conv_rate_prev = -0.1;
mean_int = zeros(N*M,N*M);
var_int = zeros(N*M,N*M);


for ite=1:n_ite
    for b=1:M*N
        mean_int_hat = zeros(taps,1);
        var_int_hat = zeros(taps,1);
        for tap_no=1:taps
            if(ind_ba(tap_no,b)~=0)
                new_chan= H(b,ind_ba(tap_no,b));
                for i2=1:1:M_mod
                    mean_int_hat(tap_no) = mean_int_hat(tap_no) + p_map(b,ind_ba(tap_no,b),i2) * alphabet(i2);
                    var_int_hat(tap_no) = var_int_hat(tap_no) + p_map(b,ind_ba(tap_no,b),i2) * abs(alphabet(i2))^2;
                end
                mean_int_hat(tap_no) = mean_int_hat(tap_no) * new_chan;
                var_int_hat(tap_no) = var_int_hat(tap_no) * abs(new_chan)^2;
                var_int_hat(tap_no) = var_int_hat(tap_no) - abs(mean_int_hat(tap_no))^2;
            end
        end
        
        mean_int_sum = sum(mean_int_hat);
        var_int_sum = sum(var_int_hat)+(no);
        for tap_no=1:taps
            if(ind_ba(tap_no,b)~=0)
                mean_int(b,ind_ba(tap_no,b)) = mean_int_sum - mean_int_hat(tap_no);
                var_int(b,ind_ba(tap_no,b)) = var_int_sum - var_int_hat(tap_no);
            end
        end
    end
    sum_prob_comp = zeros(M*N,M_mod);
    
    for a=1:M*N
        dum_sum_prob = zeros(M_mod,1);
        log_te_var=zeros(taps,M_mod);
        for tap_no=1:taps
            c=ind_ab(tap_no,a);
            if(c~=0)
                new_chan= H(c,a);               
                for i2=1:1:M_mod
                    dum_sum_prob(i2) = abs(y(c)- mean_int(c,a) - new_chan * alphabet(i2))^2;
                    dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(c,a));
                end
                dum_sum = dum_sum_prob - max(dum_sum_prob);
                dum1 = sum(exp(dum_sum));
                log_te_var(tap_no,:) = dum_sum - log(dum1);
            end
        end
        %
        % ln_qi=zeros(1,M_mod);
        for i2=1:1:M_mod
            ln_qi(i2) = sum(log_te_var(:,i2));
        end
        dum_sum = exp(ln_qi - max(ln_qi));
        dum1 = sum(dum_sum);
        sum_prob_comp(a,:) = dum_sum/dum1;
        for tap_no=1:1:taps
            c=ind_ab(tap_no,a);
            if(c~=0)
                dum_sum = log_te_var(tap_no,:);
                ln_qi_loc = ln_qi - dum_sum;
                dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                dum1 = sum(dum_sum);                
                p_map(c,a,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(c,a,:),1,M_mod);
                
            end
        end
    end
    conv_rate =  sum(max(sum_prob_comp,[],2)>0.99)/(N*M);
    if conv_rate==1
        sum_prob_fin = sum_prob_comp;
        break;
    elseif conv_rate > conv_rate_prev
        conv_rate_prev = conv_rate;
        sum_prob_fin = sum_prob_comp;
    elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
        break;
    end
     
end
for a=1:M*N
    [m1,m2]=max(sum_prob_fin(a,:));
    x_est(a)=alphabet(m2);
end
X_est=reshape(x_est,N,M).';
%% detector output
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
% scatterplot(x_data);
est_bits=reshape(qamdemod(x_data,M_mod,'gray','OutputType','bit'),N_bits_perfram,1);

end
