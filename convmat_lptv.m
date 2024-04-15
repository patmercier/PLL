%{
Title: Computationally-Efficient Linear Periodically Time-Variant 
       Digital PLL Modeling Using Conversion Matrices 
       and Uncorrelated Upsampling
Authors: Hongyu Lu, Patrick P. Mercier
Affiliation: UC San Diego
MATLAB version: R2022b update 7
System OS: Windows 11 22H2
Description: This script calls a Simulink model that simulates the output
    noise of a PLL with FPEC. Then conversion matrix algorithm is used to
    calculate the expected output phase noise and compare with the
    simulated one.
Recommand usage:
    1. Run the script with the default configuration.
    2. Change N,M,P to any desired value then run again. N will also affect
    the number of frequency domain samples = N*Npts. When using large
    N, reduce Npts accordingly.
    3. The script enables all noise sources by default. Change
    "noise_ref_en", "noise_dsm_en" and "noise_dco_en" to turn off unwanted
    noise sources.
    4. Change PLL parameters like Kp, Ki, ...
%}

% close all
%% PLL parameters
N = 16; % Reference clock Divide Ratio
M = 4; % DSM clock divid ratio
P = 16; % FPEC window length

Kp0 = 0.4; % proportional loop gain without FPEC
Kp = Kp0*N/P; % proportional loop gain
Ki = Kp0/32; % integral loop gain
Kdco = 4e6; % DCO gain [Hz/LSB]
Kpd = 300; % phase detector gain [rad^{-1}]

fref = 35e6; % reference frequency [Hz]
fdco = fref*N; % DCO output frequency [Hz]
fdsm = fref*N/M; % DSM clock frequency [Hz]

%% Configurations
noise_ref_en = 1; % 1 to enable reference noise. 0 to disable.
noise_dsm_en = 1; % 1 to enable DSM noise. 0 to disable.
noise_dco_en = 1; % 1 to enable DCO noise. 0 to disable.

Npts = 1000; % number of samples within 2pi/N. It changes how many frequence domain samples the algorithm takes
% total frequency domain samples = Npts * N
nfft_sim = 2^15; % fft length for the simulated output noise

model_name = "FPEC_sim"; % name of the Simulink model

wnorm_dt = 10e3/fdco*2*pi:2*pi/N/Npts:2*pi; % Omega, Discrete-time angular frequency
sim_time = 500*nfft_sim/fdco; % set simulation time to support 500 periodogram averages

%% Input noise sources rate conversion

tdc_white_free = 1/12/Kpd^2; % input referred tdc quantization noise power
psd_ref_up = tdc_white_free/(2*pi) /N; % upsampled TDC quantization noise

dco_white_free = 3e-5; % power of DCO input referred noise
psd_dco_up = abs(1./(1-exp(-1j*1*wnorm_dt))).^2 * dco_white_free/(2*pi) /1;

dsm_white_free = 1/12; % power of DSM quantization noise
psd_dsm_up = abs((1-exp(-1j*M*wnorm_dt)).^2).^2 * dsm_white_free/(2*pi) /M; % upsampled DSM quantization noise

%% Run Simulink simulation
% open_system(model_name); % open the simulink model
disp('Simulation starts... Takes ~30 sec')
out = sim(model_name); % run simulation

disp('Simulation finished')
disp('Conversion matrix algorithm starts...')
disp(['Number of frequency domain samples = ', num2str(length(wnorm_dt))])

%% W_N
W_NP = zeros(N,N); % Conversion matrix of the FPEC window
W_N1 = zeros(N,N); % Conversion matrix of the divider window
W_u_dsm = zeros(N,N); % Conversion matrix of DSM noise decorrelation window
W_u_ref = zeros(N,N); % Conversion matrix of REF noise decorrelation window

for i = 1:N
    for j = 1:N
        k = i-j;
        W_NP(i,j) = dtft_rect(N,P,k)/N;
        W_N1(i,j) = dtft_rect(N,1,k)/N;
        W_u_dsm(i,j) = dtft_rect_decorr(N,M,k)/N;
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
    LPF_HTF = LPF(z_vec,N,Kp,Ki,W_NP); % Conversion matrix of LPF with FPEC
    DCO_HTF = DCO(z_vec,Kdco,fdco); % Conversion matrix of DCO

    Loop_HTF = -DCO_HTF * LPF_HTF * Kpd * W_N1/N; % Conversion matrix of loop gain.

    H = (I-Loop_HTF)\I;
    H_DCO(index,:) = H(N,:); % DCO noise HTF

    H = (I-Loop_HTF)\(DCO_HTF * LPF_HTF * Kpd * W_u_ref);
    H_REF(index,:) = H(N,:); % REF noise HTF

    H = (I-Loop_HTF)\(DCO_HTF * ZOH(z_vec,M) * W_u_dsm);
    H_DSM(index,:) = H(N,:); % DSM quantization noise HTF

end

%% Calculate Output PSD
PSDtot_DCO = zeros(1,length(wnorm_dt)); % DSB output PSD
PSDtot_REF = zeros(1,length(wnorm_dt));
PSDtot_DSM = zeros(1,length(wnorm_dt));
for k = 1:N
    PSDtot_DCO = PSDtot_DCO + abs(H_DCO(:,k)').^2 .* circshift(psd_dco_up,(N-k)*Npts);
    PSDtot_REF = PSDtot_REF + abs(H_REF(:,k)').^2 .* circshift(psd_ref_up,(N-k)*Npts);
    PSDtot_DSM = PSDtot_DSM + abs(H_DSM(:,k)').^2 .* circshift(psd_dsm_up,(N-k)*Npts);
end
PSDtot_DCO_SSB = PSDtot_DCO(1:length(wnorm_dt)/2) * 2; % SSB PSD
PSDtot_REF_SSB = PSDtot_REF(1:length(wnorm_dt)/2) * 2;
PSDtot_DSM_SSB = PSDtot_DSM(1:length(wnorm_dt)/2) * 2;

PSDtot_DSB = PSDtot_DCO * noise_dco_en + PSDtot_REF * noise_ref_en + PSDtot_DSM * noise_dsm_en; % total output DSB PSD
PSDtot_SSB = PSDtot_DCO_SSB * noise_dco_en + PSDtot_REF_SSB * noise_ref_en + PSDtot_DSM_SSB * noise_dsm_en; % total output SSB PSD

%% Plot Simulated PN
[PSD_PLL,norm_f] = psd_estimate(out.PLL,nfft_sim);
semilogx(norm_f*fdco,10*log10(PSD_PLL/fdco)-3,Color=[0.4660 0.6740 0.1880])

xlim([1e4,fdco/2])
% ylim([-150,10*log10(max(PSD_PLL)/fdco)])
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

function tf = LPF(z_vec,N,Kp,Ki,W_NP)
    ZOH_N = ZOH(z_vec,N);
    P_path = W_NP*Kp*ZOH_N;
    I_path = Ki*diag(1./(1-z_vec.^(-1)));
    tf = P_path + I_path;
end

function tf = DCO(z_vec,Kdco,fdco)
    tf = diag(Kdco/fdco*z_vec.^(-1)./(1-z_vec.^(-1)));
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
