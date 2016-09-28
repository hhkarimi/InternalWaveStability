% main script to solve system of PDE describing near-inertial case
clc, clear all, close all
set(0,'defaulttextinterpreter','latex')

    % independent parameters
f = 0.1; % non-dimensional coriolis parameter
kappa = 2.36; % wavenumber of perturbation (choose maximum)
C = 0.05; alpha = f*C/2; % scaled and non-dimensionalized viscosity
sigma = 0.0; % detuning factor (sigma = sigmahat*f)

    % calculated parameters
c = sqrt(3*(1-f^2)); % group velocity
dc = 3*f; % 2nd-order dispersion factor
delta = 3*f / (2*(1-f^2)); % refraction coefficient
gamma = 3*f*sqrt(3*(1-4*f^2)) / (4*(1-f^2));

    % spatial discretization
N = 2500; eta = 50*linspace(-1,1,N); deta = (eta(end)-eta(1))/(N-1);
    
    % time discretization
dT = 0.05; % CFL: dt < deta/speed, speed set by cg and cg_eta
Tend = 400.0;
T = 0:dT:Tend;
NT = length(T);

%%%%%%%%%%%%%%%% set initial conditions, variables in (space,time)
% initial beam
Q = 1/2*exp(-eta.^2);
% initial perturbations
AMP = 10^-2;
A = AMP*Q; B(:,1) = AMP*Q;

%%%%%%%%%%%%%%%%%%% Store variables for wave plot  %%%%%%%%%%%%%%%%%%%%%%%%%
% Number of time points to plot, time steps
N2d = 10; % N2d + 1 for initial values
dtstep = round(NT / N2d);
Tplot = 0:T(dtstep):T(N2d*dtstep);
Aplot = zeros(length(eta),N2d+1); % (eta,t)
Bplot = zeros(length(eta),N2d+1);
Qplot = zeros(length(eta),N2d+1);
% initial values
Aplot(:,1) = A;
Bplot(:,1) = B;
Qplot(:,1) = Q;
iplot = 1;

% start waitbar
hw = waitbar(0,'Current Progress: 0\%');
tic
    % pass through time integration scheme
for n = 2:NT
    %%%%%%%%%%%%%%%%%%%%%  Method of lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in = [A(:); B(:); Q(:)]; % single vector input
    [tout,out] = ode45(@mlinesPDE, [T(n)-dT,T(n)], in, [], eta, sigma,c,dc,gamma,delta,alpha,kappa);
    A = out(end,1:N); B = out(end,N+1:2*N); Q = out(end,2*N+1:3*N);
    A = A(:); B = B(:); Q = Q(:); % re-shape to column vectors
        %%%% Store solution for plotting 2d waterfall
            % store values if at interval time step
    if mod(n,dtstep) == 0;
        iplot = iplot + 1;
        Aplot(:,iplot) = A;
        Bplot(:,iplot) = B;
        Qplot(:,iplot) = Q;
    end
        % update waitbar
    prog = n / NT;
    waitbar(n / NT,hw,['Current Progress: ' num2str( 100*prog ) '\%']);    
end
toc
delete(hw)
%% Plot progression as waterfall %%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaulttextinterpreter','latex')
    % decrease data points for plotting to increase plotting efficiency
eta_plot = 20*linspace(-1,1,500);
A_plot = interp1(eta,Aplot,eta_plot,'cubic',0);
B_plot = interp1(eta,Bplot,eta_plot,'cubic',0);
Q_plot = interp1(eta,Qplot,eta_plot,'cubic',0);
% eta_plot = eta; A_plot = Aplot; B_plot = Bplot; Q_plot = Qplot;

al = -26; ez = 22;

figure('name','Waterdfall 2D plot')
subplot(3,1,1)
hw = waterfall(eta_plot,Tplot,abs(Q_plot.'));
CD = get (hw, 'CData');
CD(1,:) = nan;
CD(end-2:end,:) = nan;
set(hw, 'CData', CD)
xlabel('$\eta$'); ylabel('$T$'); 
xlim([eta_plot(1) eta_plot(end)]);
zlim([0 0.5]);
% xlim([-50 50]);
title('$|Q|$'); view(al,ez);

subplot(3,1,2)
hw = waterfall(eta_plot,Tplot,abs(A_plot.'));
CD = get (hw, 'CData');
CD(1,:) = nan;
CD(end-2:end,:) = nan;
set(hw, 'CData', CD)
xlabel('$\eta$'); ylabel('$T$'); 
xlim([eta_plot(1) eta_plot(end)]);
% zlim([0 0.2]);
% xlim([-50 50]);
title('$|A|$'); view(al,ez);

subplot(3,1,3)
hw = waterfall(eta_plot,Tplot,abs(B_plot.'));
CD = get (hw, 'CData');
CD(1,:) = nan;
CD(end-2:end,:) = nan;
set(hw, 'CData', CD)
xlabel('$\eta$'); ylabel('$T$'); 
xlim([eta_plot(1) eta_plot(end)]);
% zlim([0 0.2]);
% xlim([-50 50]);
title('$|B|$'); view(al,ez);
colormap(1e-6*[1 1 1]);

%% Plot full disturbance as contour map and waterfall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % small amplitude parameter
ep = 0.1;

    % actual time
t = Tplot / ep;

    % angles and frequency
phi = asin( sigma*sqrt(ep) );
omega = 2 * sqrt( sin(phi)^2 + f^2*cos(phi)^2 );
theta = asin( sqrt( (omega^2-f^2)/(1-f^2) ));

    % coordinate orthogonal to eta
xi = linspace(eta(1),eta(end),300);
dxi = ( xi(end)-xi(1) ) / (length(xi) - 1);

    % set xy mesh
x = zeros(length(eta_plot),length(xi));
y = zeros(length(eta_plot),length(xi));
zeta = zeros(length(eta_plot),length(xi)); % wavevector direction of A and B
for n = 1:length(eta_plot)
    x(n,:) = eta_plot(n)*sin(theta) + xi.' * cos(theta);
    y(n,:) = eta_plot(n)*cos(theta) - xi.' * sin(theta);
    zeta(n,:) = x(n,:)*sin(phi) + y(n,:) * cos(phi);
end

    % pre-allocate 2D data in (eta,xi)
A2D = zeros(length(eta_plot),length(xi),length(t));
B2D = zeros(length(eta_plot),length(xi),length(t));
Q2D = zeros(length(eta_plot),length(xi),length(t));
for n = 1:length(t);
    for j = 1:length(xi)
        A2D(:,j,n) = A_plot(:,n);
        B2D(:,j,n) = B_plot(:,n);
        Q2D(:,j,n) = Q_plot(:,n);
    end
end

	% total velocity in 2D contour in x-y
% calculate stream function
    % pre-allocate
SA2D = zeros(length(eta_plot),length(xi),length(t));
SB2D = zeros(length(eta_plot),length(xi),length(t));
SQ2D = zeros(length(eta_plot),length(xi),length(t));
UA = zeros(length(eta_plot),length(xi),length(t));
UB = zeros(length(eta_plot),length(xi),length(t));
UQ = zeros(length(eta_plot),length(xi),length(t));
for n = 1:length(t) % include full wave form in eta for each time step
    SA2D(:,:,n) = ep^(3/2)/kappa * ( A2D(:,:,n) .* exp(1i*kappa*zeta/sqrt(ep) ) ) * exp(-1i * omega/2 * t(n));
    SA2D(:,:,n) = SA2D(:,:,n) + conj(SA2D(:,:,n));
    SB2D(:,:,n) = ep^(3/2)/kappa * ( B2D(:,:,n) .* exp(-1i*kappa*zeta/sqrt(ep) ) ) * exp(-1i * omega/2 * t(n));
    SB2D(:,:,n) = SB2D(:,:,n) + conj(SB2D(:,:,n));
    SQ2D(:,:,n) = ep * Q2D(:,:,n) .* exp(-1i* omega * t(n) );
    SQ2D(:,:,n) = SQ2D(:,:,n) + conj(SQ2D(:,:,n));
        % converting to along-beam velocity by centered finite differences
    U = ( SA2D(3:end,:,n) - SA2D(1:end-2,:,n) ) ./ (2*dxi);
    UA(:,:,n) = [zeros(1,length(xi)); U; zeros(1,length(xi))];
    U = ( SB2D(3:end,:,n) - SB2D(1:end-2,:,n) ) ./ (2*dxi);
    UB(:,:,n) = [zeros(1,length(xi)); U; zeros(1,length(xi))];
    U = ( SQ2D(3:end,:,n) - SQ2D(1:end-2,:,n) ) ./ (2*dxi);
    UQ(:,:,n) = [zeros(1,length(xi)); U; zeros(1,length(xi))];
end
U = UA + UB + UQ;
    % set threshhold to clean/remove lines of alternating sign on contour plots
for n = 1:length(t)
    maxc = max(max(abs(U(:,:,n))));
    dummy = U(:,:,n);
    dummy( abs(dummy) < maxc/100 ) = NaN; %maxc/1000; % can't set to 0, b/c then colors on contourf flip
    U(:,:,n) = dummy;
    maxQ = max(max(abs(UQ(:,:,n))));
    dummy = UQ(:,:,n);
    dummy( abs(dummy) < maxQ/100 ) = NaN; %maxQ/1000;
    UQ(:,:,n) = dummy;
    maxA = max(max(abs(UA(:,:,n))));
    dummy = UA(:,:,n);
    dummy( abs(dummy) < maxA/100 ) = NaN; %maxA/1000;
    UA(:,:,n) = dummy;
    maxB = max(max(abs(UB(:,:,n))));
    dummy = UB(:,:,n);
    dummy( abs(dummy) < maxB/100 ) = NaN; %maxB/1000;
    UB(:,:,n) = dummy;
end

    % prepare gray color map for contours
w_interval = 10; % even number between 0 and 100 of center
gray1 = linspace(0.8,0.6,50 - w_interval/2); % light values for negative
gray2 = ones(1,w_interval); % white space for values near 0
gray3 = linspace(0.4,0.1,50 - w_interval/2); % darker values for positive
gray = [gray1 gray2 gray3].';
gray = [gray gray gray];
    % number of contour lines
L = 7;
    % time-step at which to plot
figure('name','along-beam velocity field contour plot');
annotation('textbox', [0 0.9 1 0.1],'String',...
    ['$(kappa,\alpha,\alpha,\epsilon) =$ (',num2str(kappa),',',num2str(theta),',',num2str(alpha),...
    ',',num2str(ep),')'], ...
    'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',12);
subplot(2,2,1);
tc = 1;
[hC hC] = contourf(x,y,U(:,:,tc),L);
set(hC,'LineStyle','none');
maxc = max(max(abs(U(:,:,1))));
caxis([-maxc, maxc])
hcm = colormap(gray); hcb = colorbar;
xlabel('$x$'); ylabel('$y$');
title(['Along-beam velocity at $t/T_{\mbox{beam}} =$ ',num2str(t(tc))]);
xlim(10*[-1 1]); ylim(10*[-1 1]);

subplot(2,2,2);
tc = 3;
[hC hC] = contourf(x,y,U(:,:,tc),L);
set(hC,'LineStyle','none');
caxis([-maxc, maxc])
hcm = colormap(gray); hcb = colorbar;
xlabel('$x$'); ylabel('$y$');
title(['Along-beam velocity at $t/T_{\mbox{beam}} = $ ',num2str(t(tc))]);
xlim(10*[-1 1]); ylim(10*[-1 1]);

subplot(2,2,3);
tc = 5;
[hC hC] = contourf(x,y,U(:,:,tc),L);
set(hC,'LineStyle','none');
caxis([-maxc, maxc])
hcm = colormap(gray); hcb = colorbar;
xlabel('$x$'); ylabel('$y$');
title(['Along-beam velocity at $t/T_{\mbox{beam}} =$ ',num2str(t(tc))]);
xlim(10*[-1 1]); ylim(10*[-1 1]);

subplot(2,2,4);
tc = N2d + 1;
[hC hC] = contourf(x,y,U(:,:,tc),L);
set(hC,'LineStyle','none');
caxis([-maxc, maxc])
hcm = colormap(gray); hcb = colorbar;
xlabel('$x$'); ylabel('$y$');
title(['Along-beam velocity at $t/T_{\mbox{beam}} =$ ',num2str(t(tc))]);
xlim(10*[-1 1]); ylim(10*[-1 1]);