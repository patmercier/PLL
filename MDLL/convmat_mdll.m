%{
Title: Computationally-Efficient Linear Periodically Time-Variant 
       Digital PLL Modeling Using Conversion Matrices 
       and Uncorrelated Upsampling
Authors: Hongyu Lu, Patrick P. Mercier
Affiliation: UC San Diego
MATLAB version: R2023b
Description: This script replicate the result of [A. Santiccioli, C. Samori, 
A. L. Lacaita and S. Levantino, "Time-Variant Modeling and Analysis of 
Multiplying Delay-Locked Loops," in IEEE Transactions on Circuits and 
Systems I: Regular Papers, vol. 66, no. 10, pp. 3775-3785, Oct. 2019, 
doi: 10.1109/TCSI.2019.2918027]. 
%}

% close all
%% PLL parameters
N = 32; % Reference clock Divide Ratio

A0 = 5.7e-3; % alpha * Kpd * KT in table I

fref = 50e6; % reference frequency [Hz]
fdco = fref*N; % DCO output frequency [Hz]

%% Configurations
noise_ref_en = 1; % 1 to enable reference noise. 0 to disable.
noise_dco_en = 1; % 1 to enable DCO noise. 0 to disable.

Npts = 1000; % number of samples within 2pi/N. It changes how many frequence domain samples the algorithm takes
% total frequency domain samples = Npts * N
nfft_sim = 2^17; % fft length for the simulated output noise

model_name = "MDLL_sim"; % name of the Simulink model

wnorm_dt = 10e3/fdco*2*pi:2*pi/N/Npts:2*pi; % Omega, Discrete-time angular frequency

sim_time = 500*nfft_sim/fdco; % set simulation time to support 500 periodogram averages

%% Input noise sources rate conversion

tdc_white_free = 3e-5/10^3.73; % input referred tdc quantization noise power
% psd_ref_up = abs(1./(1-exp(-1j*N*wnorm_dt))).^2 * tdc_white_free/(2*pi)/N; % upsampled REF noise (non-white)
psd_ref_up = tdc_white_free/(2*pi) /N; % upsampled REF noise

dco_white_free = 3e-5/10^1.4; % power of DCO input referred noise
psd_dco_up = abs(1./(1-exp(-1j*1*wnorm_dt))).^2 * dco_white_free/(2*pi) /1;

%% Run Simulink simulation
% open_system(model_name); % open the simulink model
disp('Simulation starts... Takes ~60 sec')
out = sim(model_name); % run simulation

disp('Simulation finished')
disp('Conversion matrix algorithm starts...')
disp(['Number of frequency domain samples = ', num2str(length(wnorm_dt))])

%% W_N
W_N1 = zeros(N,N); % Conversion matrix of the divider window
W_u_ref = zeros(N,N); % Conversion matrix of REF noise decorrelation window

for i = 1:N
    for j = 1:N
        k = i-j;
        W_N1(i,j) = dtft_rect(N,1,k)/N;
        W_u_ref(i,j) = dtft_rect_decorr(N,N,k)/N;
    end
end

%% Transfer Functions
H_DCO = zeros(length(wnorm_dt),N); %each row contains H at different freq
H_REF = zeros(length(wnorm_dt),N);
H_DSM = zeros(length(wnorm_dt),N);

I = eye(N); % identity matrix
Omega_vec = zeros(1,N);

for index = 1:length(wnorm_dt)

    Omega_vec = wnorm_dt(index) - (N-linspace(1,N,N))*2*pi/N; % freq vector
    Omega = diag(Omega_vec);
    z_vec = exp(1j*Omega_vec); % replace z with e^{jw}
    Zm1 = diag(z_vec.^-1);
    
    A = Aofz(z_vec,A0);

    Loop_HTF = (ZOH(z_vec,N)*W_N1-I)*A*W_N1*Zm1; % Conversion matrix of loop gain.

    H = (I-Loop_HTF)\(I - ZOH(z_vec,N)*W_N1);
    H_DCO(index,:) = H(N,:); % DCO noise HTF

    H = (I-Loop_HTF)\( N*(ZOH(z_vec,N) + A - ZOH(z_vec,N)*W_N1*A) * W_u_ref); %multiply N because [20] uses time as input and output
    H_REF(index,:) = H(N,:); % REF noise HTF

end

%% Calculate Output PSD
PSDtot_DCO = zeros(1,length(wnorm_dt)); % DSB output PSD
PSDtot_REF = zeros(1,length(wnorm_dt));
% for k = N % fundamental only
% for k = 1:N-1 % aliased noise only
for k = 1:N
    PSDtot_DCO = PSDtot_DCO + abs(H_DCO(:,k)').^2 .* circshift(psd_dco_up,(N-k)*Npts);
    PSDtot_REF = PSDtot_REF + abs(H_REF(:,k)').^2 .* circshift(psd_ref_up,(N-k)*Npts);
end
PSDtot_DCO_SSB = PSDtot_DCO(1:length(wnorm_dt)/2) * 2; % SSB PSD
PSDtot_REF_SSB = PSDtot_REF(1:length(wnorm_dt)/2) * 2;

PSDtot_DSB = PSDtot_DCO * noise_dco_en + PSDtot_REF * noise_ref_en; % total output DSB PSD
PSDtot_SSB = PSDtot_DCO_SSB * noise_dco_en + PSDtot_REF_SSB * noise_ref_en; % total output SSB PSD

%% Plot Simulated PN
[PSD_PLL,norm_f] = psd_estimate(out.PLL,nfft_sim);
semilogx(norm_f*fdco,10*log10(PSD_PLL/fdco)-3,Color=[0.4660 0.6740 0.1880])

xlim([1e4,fdco/2])

hold on
grid on
box on

%% Plot Model Predicted PN
faxis = wnorm_dt(1:length(wnorm_dt)/2)/(2*pi)*fdco;
semilogx(faxis,10*log10(PSDtot_SSB*2*pi/fdco)-3,Color=[0 0.4470 0.7410])

xlabel('Frequency [Hz]')
ylabel('PN [dBc/Hz]')
legend('Simulation','Calculation')

%% Format Plot
set(gca,'linewidth',1)
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1;
end

%% Integrated Jitter
jitter_sim = std(out.PLL)/(2*pi)/fdco;
jitter_cal = sqrt(trapz(wnorm_dt,PSDtot_DSB))/(2*pi)/fdco;

disp('Conversion matrix algorithm finished')
disp(['Simulated RMS jitter = ',num2str(jitter_sim*1e12),' ps'])
disp(['Calculated RMS jitter = ',num2str(jitter_cal*1e12),' ps'])

%% Supporting Functions
function Wp0_k = dtft_rect(N,P,k)
    n = linspace(0,P-1,P);
    Wp0_k = sum(exp(-1j*n*2*pi*k/N));
end

function Wp0_k = dtft_rect_decorr(N,L,k)
    n = 0:gcd(L,N):N-1;
    Wp0_k = sqrt(gcd(L,N)) * sum(exp(-1j*n*2*pi*k/N)); 
end

function tf = ZOH(z_vec,N)
    tf = diag((1-z_vec.^(-N))./(1-z_vec.^(-1)));
end

function tf = Aofz(z_vec,A0)
    tf = diag(1./(1-z_vec.^(-1)))*diag(1./(1-z_vec.^(-1)))*A0;
end

function [PSD_SingleSide,norm_f] = psd_estimate(sig,nfft)
%   Estimate DT single-sided PSD by averaging periodograms
%   [PSD_SingleSide,norm_f] = psd_estimate(sig,nfft)
%   sig can be 1-D row or colunm vector, real valued
%   nfft is length of fft.
%   PSD_SingleSide has unit Vrms^2/1. 1 here is normalized frequency so no
%   units. When converted to CT PSD, PSD_SingleSide/fs will have unit
%   Vrms^2/Hz and norm_f*fs has unit Hz.
length_segment = nfft;
num_segments = floor(length(sig)/length_segment);
w = hann(length_segment);
U = mean(w.^2);
w_mat = repmat(w,1,num_segments);
sig_mat = reshape(sig(1:length_segment*num_segments),[length_segment,num_segments]);
Y = fft(sig_mat.*w_mat,nfft,1);
PSD = mean(abs(Y).^2/length_segment/U,2);
PSD_SingleSide = PSD(1:nfft/2+1)*2;
norm_f = linspace(0,0.5,length(PSD_SingleSide)); %normalized by fs
end
