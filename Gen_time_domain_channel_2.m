
function [G]=Gen_time_domain_channel_2(N,M,H_LS1,taps)     
G=zeros(N*M,N*M);

for q=1:N*M
    for l=1:taps
        if(mod(q,M)>=l) % due to ZP per block
            G(q,q+1-l)=H_LS1(l,q);                     % equation (42) in [R1]
        end
    end
end
end