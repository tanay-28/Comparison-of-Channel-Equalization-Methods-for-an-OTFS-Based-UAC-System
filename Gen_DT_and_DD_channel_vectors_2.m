function [nu_ml_tilda,nu_ml,K_ml]=Gen_DT_and_DD_channel_vectors_2(N,M,taps,H_LS1)
% l_max=max(L_set);
nu_ml_tilda=zeros(N,M,taps);
nu_ml=zeros(N,M,taps);
K_ml=zeros(N,N,M,taps);
Fn=dftmtx(N);
Fn=Fn./norm(Fn);
for m=1:M
    for l=1:taps
        for n=1:N
            nu_ml_tilda(n,m,l)=H_LS1(l,m+(n-1)*M);         %equation(42) in [R1]
        end
        nu_ml(:,m,l)=Fn*nu_ml_tilda(:,m,l);             %Section III-A in [R1]
        K_ml(:,:,m,l)=Fn*diag(nu_ml_tilda(:,m,l))*Fn';  %Section III-A in [R1]      
    end
end
end