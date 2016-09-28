function [R] = res_eigML(lambdahat,sigmahat,kappa,N,eta_max)
%RES Computes residuals stemming from continuous matching conditions for
%use with fsolve
%   input:
%         lambdahat  = eigenvalue
%         sigmahat  = detuning parameter
%         kappa = perturbation wavenumber
%         N  = number of nodes;
%         eta_inf = domain size;
%    output:
%         R = residual satisfying the equation of continuity
%         J = Jacobian matrix wrt to parameters l

    % run through ode solver
[~,~,~,~,errM] = odeON(lambdahat,sigmahat,kappa,N,eta_max);
% [~,~,~,~,errM] = ode_noON(l,s,k,f,N,eta_inf);
       
    % calculate residuals
R = errM;

%     % Calculate Jacobian
% D = 10^-4; % relative change to calculate Jacobian
% lJ = l*(1+D);
% [~,~,~,~,errMJ] = odeON(lJ,s,Q0,N,eta_inf);
% J = (errMJ - errM) ./ (l*D);
end

