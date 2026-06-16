function [H_t_f]=Generate_time_frequency_channel_ZP_2(N,M,H_LS1,taps)
H_t_f=zeros(N,M); % Time-frequency single tap channel matrix
Fm=dftmtx(M);
Fm=Fm./norm(Fm);
Gn=zeros(M,M);
for n=1:N
    for m=1:M
        for l=1:taps
            if(m>=l)
                Gn(m,m-l+1)=H_LS1(l,m+(n-1)*M);  %equation(42) in [R1]
            end
        end
    end
    H_t_f(n,1:M)=diag(Fm*Gn*Fm').';
end
end