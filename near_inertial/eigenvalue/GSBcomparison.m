%%%%%%% Comparison with Gerkema, Staquet, Buerot-Aubertot (2006)
clc, clear all, close all

%%% Dimensional beam characteristics from paragraphs [9] and [10]
N0 = 2*10^(-3); % buoyancy frequency (rads/sec)
U0 = 0.2; % peak velocity of beam (m/s)
omega_star = 1.405*10^(-4); % M2 forcing frequency (rads/sec)
f_star = 6.73*10^(-5); % local Coriolis parameter (rads/sec)
L_star = 1000; % estimated beam width from figure 1, (m)
% nu_star = [10^(-4),10^(-2),10^(-6)]; % viscosity: [vertical, horizontal, water], (m^2/s)
nu_star = [10^(-4),10^(-3),10^(-2)]; % viscosity: [vertical, average, horizontal], (m^2/s)

%%% nondimensional beam parameters
e = U0 / (N0*L_star);
omega = omega_star / N0;
f = f_star / N0;
theta = asin( sqrt( (omega^2-f^2) / (1-f^2) ) );
nu = nu_star ./ (N0*L_star^2);
alpha = nu ./ (2*e^2);
C = 2*alpha/f;
sigma = sqrt( ( (omega/2)^2 - f^2 ) / (e*(1-f^2)) );

%%% parameters for eigenvalue problem
sigmahat = sigma/f;

%%% results of eigenvalue problem with sigamahat above (0.95):
% lambdaBYf = [0.3873 + 0.9641i, 0.3226 + 0.9546i, 0.46 + 0.79i]; % for sigmahat = 1 and C from above
kappa = [25.67, 12.68, 2.07];
lambdaBYf = [0.81 + 0i, 0.4519, 0.3226 + 0.9546i];    

%%% spatial scale of subharmonic disturbance
k_star = kappa ./ (sqrt(e)*L_star);
wavelength = 2*pi ./ k_star;

%%% growth rate of streamfunction (velocity)
lambda_r = f*real( lambdaBYf );
growthrate = e.*lambda_r*N0; % growth rate of streamfunction perturbations 1/s
inv_gr = 1./growthrate * 1/(3600*24); % in days
inv_gr_E = inv_gr / 2; % energy growth ~ (streamfunction)^2