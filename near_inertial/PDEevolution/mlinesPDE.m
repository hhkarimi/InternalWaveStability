function [ddt] = mlinesPDE(t,in,eta,sigma,c,dc,gamma,delta,alpha,kappa)
%PDEddt returns system of ODE's in time to use in conjunction with the
% method of lines
%   sytem of ODE's in time from PDE discretized in space
%   input:
%         in = [A; B; Q];
%         A = subharmonic wave A over space (column vector)
%         B = subharmonic wave B over space (column vector)
%         Q = primary beam Q over space (column vector)
%         eta = spatial discretization
%   output:
%         ddt = right-hand side of ODEs (1 column vector comprising the 3 variables)

% organizing and discretization
N = length(in)/3;
deta = ( eta(N) - eta(1) ) / (N-1);
A = in(1:N); B = in(N+1:2*N); Q = in(2*N+1:3*N);

% derivatives (forward/backword)
% A_eta = ( 3*A(3:N) - 4*A(2:N-1) + A(1:N-2) ) ./ (2*deta);
% A_eta = [0; 0; A_eta];
% A_etaeta = ( 2*A(4:N) - 5*A(3:N-1) + 4*A(2:N-2) - A(1:N-3) ) ./ (deta^2);
% A_etaeta = [0; 0; 0; A_etaeta];
% B_eta = ( -3*B(1:N-2) - 4*B(2:N-1) + B(3:N) ) ./ (2*deta);
% B_eta = [B_eta; 0; 0];
% B_etaeta = ( 2*B(1:N-3) - 5*B(2:N-2) + 4*B(3:N-1) - B(4:N) ) ./ (deta^2);
% B_etaeta = [B_etaeta; 0; 0; 0];


%%%%%%%%%%  Centered Differenes
A_eta = ( A(3:N) - A(1:N-2) ) ./ (2*deta);
A_eta = [0; A_eta; 0];
A_etaeta = ( A(3:N) - 2*A(2:N-1) + A(1:N-2) ) ./ (deta^2);
A_etaeta = [0; A_etaeta; 0];

B_eta = ( B(3:N) - B(1:N-2) ) ./ (2*deta);
B_eta = [0; B_eta; 0];
B_etaeta = ( B(3:N) - 2*B(2:N-1) + B(1:N-2) ) ./ (deta^2);
B_etaeta = [0; B_etaeta; 0];

Q_eta = ( Q(3:N) - Q(1:N-2) ) ./ (2*deta);
Q_eta = [0; Q_eta; 0];
Q_etaeta = ( Q(1:N-2) - 2*Q(2:N-1) + Q(3:N) ) ./ (deta^2);
Q_etaeta = [0; Q_etaeta; 0];


dAdt = -sigma*c/kappa * A_eta + 1i*dc/(2*kappa^2) * A_etaeta - 2*alpha*kappa^2 * A ...
    + 1i*delta*kappa^2*abs(Q_eta).^2 .* A - gamma*Q_etaeta .* conj(B);
dBdt = sigma*c/kappa * B_eta + 1i*dc/(2*kappa^2) * B_etaeta - 2*alpha*kappa^2 * B ...
    + 1i*delta*kappa^2*abs(Q_eta).^2 .* B - gamma*Q_etaeta .* conj(A);
dQdt = -2*gamma * (A.*B);
% dQdt = zeros(N,1); % linearized case

ddt = [dAdt; dBdt; dQdt];
end