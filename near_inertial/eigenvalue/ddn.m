function [dydeta] = ddn(eta,y,lambdahat,sigmahat,kappa,deta)
%DDN system of differential equations
%   input:
%         eta = evaluation point
%         y   = initial conditions
%         lambdahat   = eigenvalue (parameter)
%         sigmahat   = detuning parameter
%         kappa  = pertubation wavenumber
%         deta = spatial discretization (uniform mesh, required for varying coefficients)
%   output:
%         dydeta = [dAdeta; d2Adeta2; dBdeta; d2Bdeta2]

    % primary beam
Qb = beam_profile(eta-deta);
Qc = beam_profile(eta); % beam profile
Qf = beam_profile(eta+deta);

        % derivatives
Q_eta = ( Qf-Qb ) ./ (2*deta);
Q_etaeta = ( Qf - 2*Qc + Qb ) ./ (deta^2);

% %%%%%%%%%%%%%%%   Test if the beam profile is correctly loading
% figure(10)
% subplot(3,1,1); plot(eta,abs(Fc),'.'); hold on;
% subplot(3,1,2); plot(eta,abs(Q_eta/Q0),'.'); hold on;
% subplot(3,1,3); plot(eta,abs(Q_etaeta/Q0),'.'); hold on;

        % state-space ODE
dydeta = [y(2); % Ahat_eta
          -2i*kappa^2/3*lambdahat*y(1) - 2i*sigmahat*kappa/sqrt(3)*y(2) ...
            - kappa^4*abs(Q_eta).^2*y(1) - 1i*sqrt(3)*kappa^2/2*Q_etaeta*y(3); % Ahat_etaeta
          y(4); % Bhat_eta_star
          2i*kappa^2/3*lambdahat*y(3) - 2i*sigmahat*kappa/sqrt(3)*y(4) ...
            - kappa^4*abs(Q_eta).^2*y(3) + 1i*sqrt(3)*kappa^2/2*Q_etaeta*y(1)]; % Bhat_etaeta_star
end