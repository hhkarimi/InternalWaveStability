%%%%  Find lambdahat for fixed sigmahat and sweeps kappa up or down
clc, clear all, close all
    % set default interpreter to latex
set(0,'defaulttextinterpreter','latex')

    % load parameters (k,sf) and initial guesses (dnew)
%[filename, pathname] =  uigetfile('*.mat','Load parameters and initial guesses for lambdahat');
%completename = fullfile(pathname, filename);
%load(completename);
    % loads lambdahat, sigmahat, and kappa: ALSO CHANGE IN SAVED FILE NAME WHEN CHOOSING A DIFFERENT EV BRANCH!
% lambdahat1 = lambdahat(1); % ALSO CHANGE IN SAVED FILE NAME WHEN CHOOSING A DIFFERENT EV BRANCH!
% lambdahat1 = lambdahat(2); % ALSO CHANGE IN SAVED FILE NAME WHEN CHOOSING A DIFFERENT EV BRANCH!
lambdahat1 = 0.15;
kappa_ind = linspace(kappa,30,1000);

    % set spatial discretization and  initial guesses
N = 5000;
eta_max = 5;
eta = eta_max*linspace(-1,1,N);
deta = (eta(end) - eta(1)) / (N-1);
        
maxit = 30; % max iterations
es = 10^-3; % percent error
    % pre-allocate
lambdahat = NaN(1,length(kappa_ind));
    % Options with fsolve
options=optimset('MaxIter',1e3,'TolFun',1e-3);
hw = waitbar(0,'Current Progress: 0\%');
    % initialize flags
iter = 0; exitflag = 1; 
tic
for i = 1:length(kappa_ind)    
        % initialize first set of guesses
    if i == 1; lambdahat(i) = lambdahat1;
    else lambdahat(i) = lambdahat(i-1);
    end
    
        % solve for lambdahat if solutions may exist
    if isfinite(lambdahat(i)) == 1;
                % with Newton-Raphson
%         [l1(i),~,ea,iter] = newtmult(@res_eig,l1(i),es,maxit,s,Q0_ind(i),N,eta_max);
                % with fsolve
        [lambdahat(i),~,exitflag,exitout] = fsolve(@(l) res_eigML(l,sigmahat,kappa_ind(i),N,eta_max),lambdahat(i),options);
    else lambdahat(i) = NaN*(1+1i);
    end
    if iter == maxit; keyboard; lambdahat(i) = NaN*(1+1i); iter = 0; end
    if exitflag ~= 1; keyboard; lambdahat(i) = NaN*(1+1i); exitflag = 1; end
        % update waitbar
    prog = i / length(kappa_ind);
    waitbar(i / length(kappa_ind),hw,['Current Progress: ' num2str( 100*prog ) '\%']);
end
toc
delete(hw)

%%    % Plot eigenvalue solutions
set(0,'defaulttextinterpreter','latex')
figure(1)
plot(kappa_ind,real(lambdahat)); hold on; plot(kappa_ind,imag(lambdahat),'--r');
grid on
title(['Eigenvalues for $\hat{\lambda}(\hat{\sigma},\kappa), \quad \hat{\sigma} =$ ',...
    num2str(sigmahat)]);
xlabel('$k$','interpreter','latex'); 
ylabel('$\Lambda$','interpreter','latex');
legend('real','imag')
%% saving data and keyboard
if kappa > kappa_ind(end); direction = 'down'; else direction = 'up'; end
savefile = sprintf('sh%dkappa%.1flh2%s.mat',sigmahat,kappa,direction);
save(savefile,'sigmahat','kappa_ind','lambdahat');
break

%% generate mode shapes for some kappa as diagnostic
set(0,'defaulttextinterpreter','latex')
    % initial values for mode shape at chosen kappa
% ind_check = 2;
% l0 = lambdahat(ind_check);
% kappa_check = kappa_ind(ind_check);
sigmahat = 0;
l0 = 0.3202 + 0.7693i; kappa_check = 2.364;
N = 4000; eta_max = 20;

% Solve ODE with refined parameter
[ym,etam,yp,etap,errM] = odeON(l0,sigmahat,kappa_check,N,eta_max);

    % normalize eigenfunctions by setting = 1 at eta ~ 0
Y = ym(end,1);
ym(:,1) = ym(:,1) / Y;
yp(:,1) = yp(:,1) / Y;
ym(:,3) = ym(:,3) / Y;
yp(:,3) = yp(:,3) / Y;
figure(2);
subplot(2,2,1)
plot(etam,real(ym(:,1))); hold on; plot(etam,imag(ym(:,1)),'--r');
plot(etap,real(yp(:,1))); plot(etap,imag(yp(:,1)),'--r');
xlabel('$\eta$'); hl = legend('real $\hat{A}$','imag $\hat{A}$');
set(hl,'interpreter','latex')
title(['$\hat{\lambda} =$ ',num2str(l0),' at $k$ =' ,num2str(kappa_check)]);
subplot(2,2,2)
plot(etam,real(ym(:,3))); hold on; plot(etam,imag(ym(:,3)),'--r');
plot(etap,real(yp(:,3))); plot(etap,imag(yp(:,3)),'--r');
xlabel('$\eta$'); hl = legend('real $\hat{B}$','imag $\hat{B}$');
set(hl,'interpreter','latex')
subplot(2,2,[3 4]); 
plot(etam,abs(ym(:,1))); hold on; plot(etap,abs(yp(:,1)));
plot(etam,abs(ym(:,3)),'--r'); hold on; plot(etap,abs(yp(:,3)),'--r');
