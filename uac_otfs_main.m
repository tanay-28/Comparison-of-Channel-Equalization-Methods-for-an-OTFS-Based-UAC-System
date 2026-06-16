close all
clear all

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 16;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
SNR_dB = 0:5:40;
% SNR_dB=30;
errors_MRC=zeros(1,length(SNR_dB));
errors_MPA=zeros(1,length(SNR_dB));
errors_LMMSE=zeros(1,length(SNR_dB));
errors_1tap=zeros(1,length(SNR_dB));
%% delay-Doppler grid symbol placement
% ZP length  should be set to greater than or equal to maximum value
% of delay_taps
length_ZP = M/16;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-length_ZP;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;



% Time and frequency resources
B=6400;
df=B/M;
dt=0.1;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
%% plot channel key
plotting_key =0;
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
Fm=dftmtx(M);
Fm=Fm./norm(Fm);
% current_frame_number=zeros(1,length(SNR_dB));
T = 1;
for t=1:T
for iesn0 = 1:length(SNR_dB)
    
        
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        X = Generate_2D_data_grid(N,M,data,data_grid);
        
        
        %% OTFS modulation%%%%
        X_tilda=X*Fn';
        s = reshape(X_tilda,N*M,1);
        
        
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code. 
        %% synthetic channel model with equal power paths with delays [0,l_max] and Dopplers [-k_max,k_max]
%         taps=4;
%         l_max=delay_spread;
%         k_max=4;
%         chan_coef=1/sqrt(2)*(randn(1,taps)+1i.*randn(1,taps));
%         delay_taps=randi(taps,[1,l_max+1]);  
%         delay_taps=sort(delay_taps-min(delay_taps));  %% random delay shifts in the range [0,l_max] 
%         Doppler_taps=k_max-2*k_max*rand(1,taps);   %% uniform Doppler profile [-k_max,k_max]
%         L_set=unique(delay_taps);

        
              %%
              % close all;

%% Parameters to change
h0=100; % surface height (depth) [m]
ht0=40; % TX height [m]
hr0=20; % RX height [m]
d0=1000; % channel distance [m]
%%
%5250
k=1.7; % spreading factor
c=1500; % speed of sound in water [m/s]
c2=1200; % speed of sound in bottom [m/s] (>1500 for hard, < 1500 for soft)
cut=20; % do not consider arrivals whose strength is below that of direct arrival divided by cut

fmin=10e3; % minimum frequency [Hz]
% B=6400; % bandwidth [Hz]
%% Small-Scale (S-S) parameters: 
sig2s= 1.125; % variance of S-S surface variations 
sig2b=sig2s/2; % variance of S-S bottom variations 
 
B_delp= 5e-4; % 3-dB width of the p.s.d. of intra-path delays (assumed constant for all paths)
Sp= 20; % number of intra-paths (assumed constant for all paths)
mu_p= .5/Sp; % mean of intra-path amplitudes (assumed constant for all paths)
nu_p= 1e-6; % variance of intra-path amplitudes (assumed constant for all paths)

%% Large-Scale (L-S) parameters:
T_tot=204.8;%seconds
t_tot_vec=(0:dt:T_tot-dt); Lt_tot=length(t_tot_vec);

h_bnd=[-10 10]; % range of surface height variation (L-S realizations are limited to h+h_band)
ht_bnd=[-5 5]; % range of transmitter height variation
hr_bnd=[-5 5]; % range of receiver height variation
d_bnd=[-20 20]; % range of channel distance variation
sig_h=1; % standard deviation of L-S variations of surface height 
sig_ht=1; % standard deviation of L-S variations of transmitter height 
sig_hr=1; % standard deviation of L-S variations of receiver height 
sig_d=1; % standard deviation of L-S variations of distance height 
a_AR= .9; % AR parameter for generating L-S variations (constant for variables h, ht, hr, d) 

method='sim_seq';

%% Doppler: 
% All parameters can be entered as scalars or vectors of size(t_tot_vec).

% drifting: 
vtd= 0.1*cos((1/20)*(2*pi*t_tot_vec)-rand*2*pi); %transmitter drift velocity,what is 1/20?
theta_td= linspace(rand*2*pi,rand*2*pi,length(t_tot_vec)); 

vrd= 0.02*cos((1/5)*(2*pi*t_tot_vec)); %receiver drift velocity,what is 1/5?
theta_rd= rand*2*pi; 

% vehicular: 
vtv= interp1(0:T_tot/4:T_tot, 1*randn(size(0:T_tot/4:T_tot)), t_tot_vec, 'spline');   
% vtv=0;
theta_tv= 0; 
vrv= 0;  
theta_rv= 0; 

% surface: 
Aw= .05; fw= .01; 


Dopp_mat= zeros(Lt_tot, 10);
if length(vtd)==1; Dopp_mat(:,1)= repmat(vtd, Lt_tot, 1); 
else Dopp_mat(:,1)= vtd.';
end
if length(theta_td)==1; Dopp_mat(:,2)= repmat(theta_td, Lt_tot, 1); 
else Dopp_mat(:,2)= theta_td.';
end
if length(vrd)==1; Dopp_mat(:,3)= repmat(vrd, Lt_tot, 1); 
else Dopp_mat(:,3)= vrd.';
end
if length(theta_rd)==1; Dopp_mat(:,4)= repmat(theta_rd, Lt_tot, 1); 
else Dopp_mat(:,4)= theta_rd.';
end
if length(vtv)==1; Dopp_mat(:,5)= repmat(vtv, Lt_tot, 1); 
else Dopp_mat(:,5)= vtv.';
end
if length(theta_tv)==1; Dopp_mat(:,6)= repmat(theta_tv, Lt_tot, 1); 
else Dopp_mat(:,6)= theta_tv.';
end
if length(vrv)==1; Dopp_mat(:,7)= repmat(vrv, Lt_tot, 1); 
else Dopp_mat(:,7)= vrv.';
end
if length(theta_rv)==1; Dopp_mat(:,8)= repmat(theta_rv, Lt_tot, 1); 
else Dopp_mat(:,8)= theta_rv.';
end
if length(Aw)==1; Dopp_mat(:,9)= repmat(Aw, Lt_tot, 1); 
else Dopp_mat(:,9)= Aw.';
end
if length(fw)==1; Dopp_mat(:,10)= repmat(fw, Lt_tot, 1); 
else Dopp_mat(:,10)= fw.';
end


%%
fmax=fmin+B; 
f_vec=(fmin:df:fmax-df).'; Lf=length(f_vec);
fc=(fmax+fmin)/2; 
f0=fmin;
f_vec2=(f0-B/2:df:f0+B+B/2); dif_f=f_vec2-fc;
T_SS=T_tot/4;
t_vec=(0:dt:T_SS-dt); Lt=length(t_vec);

%% Large-scale channel parameters:
N_LS=round(T_tot/T_SS); 
t_tot_vec=(0:dt:T_tot-dt); Lt_tot=length(t_tot_vec);
 

%% Large-Scale (L-S) & Small-Scale (S-S) Simulation Methods:
method_LS= lower(method(1:3));
method_SS= lower(method(5:7));

%% Doppler parameters: 

% % drifting: 
vtd_tot= Dopp_mat(:,1).'; 
theta_td_tot= Dopp_mat(:,2).'; 
vrd_tot= Dopp_mat(:,3).';
theta_rd_tot= Dopp_mat(:,4).';

% vehicular: 
vtv_tot= Dopp_mat(:,5).';
theta_tv_tot= Dopp_mat(:,6).';
vrv_tot= Dopp_mat(:,7).';
theta_rv_tot= Dopp_mat(:,8).';

% surface: 
Aw_tot= Dopp_mat(:,9).';
fw_tot= Dopp_mat(:,10).'; 

%% Large-scale loop: 
H_LS= zeros(Lf, Lt*N_LS);
del_h=0; del_ht=0; del_hr=0; del_d=0; 
h=h0; ht=ht0; hr=hr0; d= d0; % initial values

%% initialize plotting of large-scale channel parameters: 
adopp0=zeros(1,50);

for LScount= 1:N_LS
rndvec= randn(1, 4); 


del_h= a_AR * del_h + sqrt(1-a_AR^2)*sig_h*rndvec(1); 
if (del_h > h_bnd(2))||(del_h < h_bnd(1)), 
    del_h= del_h- 2*sqrt(1-a_AR^2)*sig_h*rndvec(1); 
end
htemp=h;
h= h0+del_h;

del_ht= a_AR * del_ht + sqrt(1-a_AR^2)*sig_ht*rndvec(2); 
if (del_ht > ht_bnd(2))||(del_ht < ht_bnd(1)), 
    del_ht= del_ht- 2*sqrt(1-a_AR^2)*sig_ht*rndvec(2); 
end
httemp=ht;
ht= ht0+del_ht;


del_hr= a_AR * del_hr + sqrt(1-a_AR^2)*sig_hr*rndvec(3); 
if (del_hr > hr_bnd(2))||(del_hr < hr_bnd(1)), 
    del_hr= del_hr- 2*sqrt(1-a_AR^2)*sig_hr*rndvec(3); 
end
hrtemp=hr;
hr= hr0+del_hr;

del_d= a_AR * del_d + sqrt(1-a_AR^2)*sig_d*rndvec(4); 
if (del_d > d_bnd(2))||(del_d < d_bnd(1)), 
    del_d= del_d- 2*sqrt(1-a_AR^2)*sig_d*rndvec(4); 
end
dtemp=d;
d= d0+del_d;
end

% Find Large-scale model parameters: 

switch method_LS
case 'sim'
    [lmean,taumean,Gamma,theta,ns,nb,hp]= mpgeometry(h,h-ht,h-hr,d,fc,k,cut,c,c2);

case 'bel'
    [lmean,taumean,theta,ns,nb,hp,p0_ind]=runBellhop(h,h-ht,h-hr,d,fc,c,c2);
    lmean=[lmean(p0_ind), lmean(1:p0_ind-1), lmean(p0_ind+1:end)];
    taumean=[taumean(p0_ind), taumean(1:p0_ind-1), taumean(p0_ind+1:end)];
    theta=[theta(p0_ind), theta(1:p0_ind-1), theta(p0_ind+1:end)];
    ns=[ns(p0_ind), ns(1:p0_ind-1), ns(p0_ind+1:end)];
    nb=[nb(p0_ind), nb(1:p0_ind-1), nb(p0_ind+1:end)];
    hp=[hp(p0_ind), hp(1:p0_ind-1), hp(p0_ind+1:end)];

end

lmean= lmean(taumean<1/df); 
 theta= theta(taumean<1/df); 
 ns= ns(taumean<1/df); 
  nb= nb(taumean<1/df); 
  hp= hp(taumean<1/df); 
  taumean= taumean(taumean<1/df);
one_delay_tap=1/B;
taumean1=round(taumean/one_delay_tap);
L_set=unique(taumean1);
  l_max=max(taumean1);
  taps=length(taumean);
% P=l_max;
H1 = zeros(taps,Lt*Lf);
H_LS1 = zeros(taps,Lf*Lt_tot);
for LScount= 1:N_LS
rndvec= randn(1, 4); 
% ignore paths with delays longer than allowed by frequency resolution:  

 P=length(taps); % number of paths



% Reference path transfer function: 
H0= 1./sqrt(lmean(1)^k* (10.^(absorption(f_vec/1000)/10000)).^lmean(1) );
H2= hp(1)*repmat( exp(-1j*2*pi*f_vec*taumean(1)) , 1, Lt);
% H_int=hp(1)*repmat( exp(-1j*2*pi*f_vec*taumean(1)) , 1, Lt);
%% Find small-scale model parameters:

sig_delp= sqrt(1/c^2*((2*sin(theta)).^2).*(ns*sig2s+nb*sig2b)); 
rho_p= exp(-((2*pi*f_vec).^2) * (sig_delp.^2/2));
rho_p_bb= exp(-(2*pi*(f_vec-fmin).^2) * (sig_delp.^2/2));


Bp= ((2*pi*f_vec*sig_delp).^2).*B_delp;
sig_p= sqrt(.5*(mu_p.^2.*Sp.*(1-rho_p.^2)+Sp.*nu_p.^2));

%% Find doppler rates:
% quarter section of orginal vectors of size 256 ie. 64
% drifting: 
vtd= vtd_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1)); 
theta_td= theta_td_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
vrd= vrd_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
theta_rd= theta_rd_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));

% vehicular: 
vtv= vtv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
theta_tv= theta_tv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
vrv= vrv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
theta_rv= theta_rv_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
 
% surface: 
Aw= Aw_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
fw= fw_tot(1+(LScount-1)*(Lt-1):1+LScount*(Lt-1));
vw= 2*pi*fw.*Aw;


% first path doppler:

vdrift= vtd.*cos(theta(1)-theta_td)-vrd.*cos(theta(1)+theta_rd); 
adrift= vdrift/c; 

vvhcl= 0; 
avhcl= vvhcl/c;

vsurf=  0; 
asurf= vsurf/c; 


adopp= adrift + avhcl + asurf*ns(1); 
eff_adopp = adopp0(1)+cumsum(adopp);
Dopp= exp(1j*2*pi*f_vec*(eff_adopp*dt)); 
adopp0(1)=eff_adopp(end);
H2= H2.*Dopp;
% H_int=H_int.*Dopp;

%% small-scale simulation:

for p=2:length(lmean) 

switch method_SS
    case 'dir' % DIRectly generate gamma: 

    gamma=zeros(Lf, Lt);
    for counti=1:Sp
        % randsum=0;
        % randsum1=0;
        % for t=1:100
        %     randsum=randsum+randn(1,Lt);
        %     randsum1=randsum1+randn(1,2*Lt);
        % end
        % randsum=randsum/100;
        % randsum1=randsum1/100;
        gamma_pi= mu_p+nu_p*randn(1,Lt); 
        gamma_pi=repmat(gamma_pi, Lf, 1)*Sp; 
        deltau_pi=zeros(Lf, Lt);
        w_delpi= sig_delp(p)*sqrt(1-exp(-1*pi*B_delp*dt)^2)*randn(1,2*Lt);
        temp_deltau_pi=filter(1, [1, -exp(-pi*B_delp*dt)], w_delpi);
        for countf= 1:Lf
            deltau_pi(countf, :)= temp_deltau_pi(Lt+1:end);
        end
        gamma=gamma+ gamma_pi.*exp(-1j*2*pi*repmat(f_vec, 1,Lt).*deltau_pi);
    end
    
case 'seq' % use the Statistically EQuivalent model:
    if sig_delp(p)*f0 < .4
    fprintf(['\n Warning! sig_delp*f0 < 0.4 for path ' num2str(p),' of L-S realization ', num2str(LScount),'. The statistically equivalent model is not accurate. \n'])
    end
    alpha_p= exp(-pi*Bp(:,p)*dt); 
    Alpha_p= diag(alpha_p);    
    
    Wp_cov= 2*(ones(Lf, Lf)-alpha_p*alpha_p.').*toeplitz(rho_p_bb(:,p)).*(sig_p(:,p)*sig_p(:,p).');

    error_takagi=-1;
    while sign(error_takagi)~=1
        [Wp_Q, Wp_s, error_takagi]= find_takagi_factor(Wp_cov);  
    end
    
    Wp_factor= Wp_Q*sqrt(diag(Wp_s)); 
    wp= (Wp_factor)*(randn(Lf,Lt) + 1j*randn(Lf,Lt));
        
    gammabar=repmat(mu_p + mu_p*Sp*rho_p(:,p), 1, Lt);

    Delgamma=zeros(Lf,Lt+1);
    for count_t=1:Lt
        Delgamma(:,count_t+1)= Alpha_p*Delgamma(:,count_t)+wp(:,count_t);
    end
    Delgamma=Delgamma(:, 2:end);

   gamma= gammabar+Delgamma; 

end

    %% Doppler term:
    vdrift= vtd.*cos(theta(p)-theta_td)-vrd.*cos(theta(p)+theta_rd); 
    adrift= vdrift/c; 
    
    vvhcl= vtv.*cos(theta(p)-theta_tv)-vrv.*cos(theta(p)+theta_rv)-(vtv.*cos(theta(1)-theta_tv)-vrv.*cos(theta(1)+theta_rv));
    avhcl= vvhcl/c;

    phi_pj= 2*pi*rand(1,ns(p))-pi;
    sum_j= zeros(1, Lt); 
    for jcount= 1: ns(p)
        sum_j= sum_j + sin(phi_pj(jcount)+2*pi*fw.*t_vec);
    end
    vsurf= 2*vw.*sin(theta(p)).*sum_j;
    asurf= vsurf/c; 
    
    adopp= adrift + avhcl + asurf*ns(p);
    eff_adopp = adopp0(p)+cumsum(adopp);
    Dopp= exp(1j*2*pi*f_vec*eff_adopp*dt); 
    adopp0(p)=eff_adopp(end);
    % If adopp was constant, we would have: 
    % Dopp= exp(1j*2*pi*repmat(adopp, Lf, 1).*repmat(f0, Lf, Lt));     
    
    % Multiply gamma by hp:
    gamma_noDopp=gamma;
    gamma=gamma.*Dopp;
    H2= H2+ hp(p)*repmat( exp(-1j*2*pi*f_vec*taumean(p)) , 1, Lt).*gamma;
    H_int = hp(p)*repmat( exp(-1j*2*pi*f_vec*taumean1(p)) , 1, Lt).*gamma;
    H1(p,:) = reshape(H_int,[1, Lt*Lf]);

end

H2= repmat(H0 ,1, Lt).*H2;
H_LS(:, (LScount-1)*Lt+1:(LScount)*Lt)= H2; 
H_LS1(:, ((LScount-1)*Lt)*Lf+1:(LScount)*Lf*Lt)= H1; 
H_LS1 = ifft(H_LS1);
end
Lt_tot= size(H_LS1, 2);
hmat=zeros(Lf,Lt*4);
for countt=1:Lt*4
    hmat(:, countt)= ifft(H_LS(:, countt)); %hmat(:, countt)=fftshift(hmat(:, countt));
end
% for countt=1:Lt*4
% hmat(:,countt)=(Fm)*(H_LS(:,countt));
% end
% hmat=fftshift(hmat*(dftmtx(M*N)^-1));
if plotting_key ==1
shift=10; skip=10; 
figure; axes('fontsize', 16);
hmat2plot= abs(hmat(1:end, 1:skip:end)).';
image(((0:Lf-1)-shift)/B*1000, (0:skip:Lt_tot-1)*dt, (circshift(abs(hmat(1:end, 1:skip:end)), shift)).',  'CDataMapping','scaled');
xlabel('delay [ms]', 'fontsize', 16), ylabel('time [s]', 'fontsize', 16), axis('ij'); 
set(gca,'YTick', 0:T_SS:T_tot);
colorbar;
end
        %% channel output%%%%%
        [G]=Gen_time_domain_channel_2(N,M,H_LS1,taps);
        [H,H_tilda,P]= Gen_DD_and_DT_channel_matrices(N,M,G,Fn);
        r=zeros(N*M,1);
        noise= sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
        % l_max=max(delay_taps);
        for q=1:N*M
            for l=1:taps
                if(q>=l)
                    r(q)=r(q)+H_LS1(l,q)*s(q-l+1);
                end
            end
        end
        r=r+noise;
        
        %% OTFS demodulation%%%%
        Y_tilda=reshape(r,M,N);
        Y = Y_tilda*Fn;
        y=reshape(Y.',N*M,1);
       

        %% test: the received time domain signal can be generated element in the matrix form (using r=G.s+noise).
        %         r_test=G*s+noise;
        %         test_delay_time_matrix_error=norm(r_test-r)
        %% test: the received DD domain signal can be generated in the matrix form (using y=H.x+noise).
        %         noise_DD=kron(eye(M),Fn)*P'*noise;
        %         x=reshape(X.',N*M,1);
        %         y_test=H*x_vec+noise_DD;
        %         text_delay_Doppler_matrix_error=norm(y_test-y)
        
        %% Generate delay-time and delay-Doppler channel vectors from the time domain channel.
%         [nu_ml_tilda]=Gen_delay_time_channel_vectors(N,M,l_max,gs);
         [nu_ml_tilda,nu_ml,K_ml]=Gen_DT_and_DD_channel_vectors_2(N,M,taps,H_LS1); 
        
        %% Generate block-wise time-frequency domain channel
         [H_tf]=Generate_time_frequency_channel_ZP_2(N,M,H_LS1,taps);
        
        %% Different detection methods
        
        n_ite_MRC=15; % maximum number of MRC detector iterations  (Algorithm 2 in [R1])
        n_ite_algo3=15; % maximum number of matched filtered Guass Seidel (Algorithm 3 in [R1])
        n_ite_MPA=15; % maximum number of MPA detector iterations
        %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        omega=1.25;
        if(M_mod>=64)
            omega=0.25;     % set omega to a smaller value (for example: 0.05) for modulation orders greater than 64-QAM
        end
        decision=1; %1-hard decision, 0-soft decision
        init_estimate=1; %1-use the TF single tap estimate as the initial estimate for MRC detection, 0-initialize the symbol estimates to 0 at the start of MRC iteration
        %(Note: it is recommended to set init_estimate to 0 for higher order modulation schemes like 64-QAM as the single tap equalizer estimate may be less accurate)
        
        %MRC detectors in [R1]
        
        [est_info_bits_MRC,det_iters_MRC,data_MRC] = MRC_delay_time_detector_21(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,nu_ml_tilda,taps,omega,decision,init_estimate,n_ite_MRC);
         [est_info_bits_MPA,det_iters_MPA,data_MPA] = MPA_detector_21(N,M,M_mod,sigma_2(iesn0),data_grid,y,H,n_ite_MPA);
         [est_info_bits_1tap,data_1tap] = TF_single_tap_equalizer_21(N,M,M_mod,sigma_2(iesn0),data_grid,Y,H_tf);
         [est_info_bits_LMMSE,data_LMMSE] = Block_LMMSE_detector_21(N,M,M_mod,sigma_2(iesn0),data_grid,r,H_LS1,taps);
        
        
        %% errors count%%%%%
        errors_MRC(1,iesn0) = errors_MRC(1,iesn0)+sum(xor(est_info_bits_MRC,trans_info_bit))/N_bits_perfram;
        errors_MPA(1,iesn0) = errors_MPA(1,iesn0)+ sum(xor(est_info_bits_MPA,trans_info_bit))/N_bits_perfram;
        errors_1tap(1,iesn0) = errors_1tap(1,iesn0)+sum(xor(est_info_bits_1tap,trans_info_bit))/N_bits_perfram;
        errors_LMMSE(1,iesn0) = errors_LMMSE(1,iesn0)+sum(xor(est_info_bits_LMMSE,trans_info_bit))/N_bits_perfram;
        
    end
    end   
  
    

errors_MRC(1,iesn0)=errors_MRC(1,iesn0)./T;
errors_MPA(1,iesn0)=errors_MPA(1,iesn0)./T;
errors_1tap(1,iesn0)= errors_1tap(1,iesn0)./T;
errors_LMMSE(1,iesn0)=errors_LMMSE(1,iesn0)./T;
figure;
% semilogy(SNR_dB,avg_ber_Algo1,'-x','LineWidth',2,'MarkerSize',8)
% hold on
semilogy(SNR_dB,errors_MRC,'-s','LineWidth',2,'MarkerSize',8)
hold on
% semilogy(SNR_dB,avg_ber_Algo3,'-x','LineWidth',2,'MarkerSize',8)
% hold on
% semilogy(SNR_dB,avg_ber_Algo3_low_complexity,'-x','LineWidth',2,'MarkerSize',8)
% hold on
semilogy(SNR_dB,errors_MPA,'-s','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,errors_1tap,'-s','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,errors_LMMSE,'-s','LineWidth',2,'MarkerSize',8)
legend('MRC','MPA','time-freq single tap','time-domain LMMSE')
grid on
xlabel('SNR(dB)')
ylabel('BER')

