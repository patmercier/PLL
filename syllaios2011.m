%{
Title: Linear Periodically Time-Variant Digital PLL Phase
       Noise Modeling Using Conversion Matrices and
       Uncorrelated Upsampling
Authors: Hongyu Lu, Patrick P. Mercier
Affiliation: UC San Diego
MATLAB version: R2023b
Description: This script replicate the result of [I. L. Syllaios and 
P. T. Balsara, "Multi-clock domain analysis and modeling of all-digital 
frequency synthesizers," 2011 IEEE International Symposium of Circuits 
and Systems (ISCAS), Rio de Janeiro, Brazil, 2011, pp. 153-156, 
doi: 10.1109/ISCAS.2011.5937524.]. Only DSM noise is implemented. 

Before using this script, run convmat_fpec.m first to set up the variables.
make the following changes in convmat_fpec.m:
1. Set P = N to disable FPEC.
2. Set noise_ref_en and noise_dco_en to 0.
%}

H_DSM = zeros(length(wnorm_dt),N);

for index = 1:length(wnorm_dt)
Omega_vec = wnorm_dt(index) - (N-linspace(1,N,N))*2*pi/N; % freq vector
z = exp(1j*Omega_vec);
HD = (1-z.^(-M)).^2;
HR1 = (1-z.^(-N))./(1-z.^(-1));
HR2 = (1-z.^(-M))./(1-z.^(-1));
F = (Kp0 + Ki * 1./(1-z.^(-N)))*Kpd/N;
HI = Kdco/fdco./(z-1);

Halias = F.*HR1.*HI./(N+F.*HR1.*HI);
Hndsm = N.*HD.*HR2.*HI./(N+F.*HR1.*HI);
A = 1./(1-Halias);

H_DSM(index,:) = abs(Halias(N))^2 * abs(Hndsm).^2 .* abs(A).^2;
H_DSM(index,N) = abs(Hndsm(N))^2;
end

PSDtot_DSM = zeros(1,length(wnorm_dt));
for k = 1:N
    PSDtot_DSM = PSDtot_DSM + H_DSM(:,k)' .* dsm_white_free/(2*pi) /M;
end

PSDtot_DSM_SSB = PSDtot_DSM(1:length(wnorm_dt)/2) * 2;
semilogx(faxis,10*log10(PSDtot_DSM_SSB*2*pi/fdco)-3,LineWidth=2)
hold on
