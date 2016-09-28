function [ym,etam,yp,etap,errM] = odeON(lambdahat,sigmahat,kappa,N,eta_max)
%ODEON uses ode45 to solve ODE with orthonormalization
%     input:
%        lambdahat = eigenvalue
%        sigmahat = detuning parameter
%        N = number of nodes;
%        eta_max = domain size;
%     output:
%        ym = solution at etam marched from left
%        yp = solution at etap marched from right

    % orthonormalization criteria
a_max = 8/9*pi; % angle between solution spaces

    % set spatial coordinates
eta = linspace(-eta_max,eta_max,N); % must be uniform for convenient in computing derivatives of beam profile
deta = eta(2)-eta(1);

    % march solution for different initial values starting far from beam (asymptotic, decoupled solutions)
alpha(1) = -1i*sigmahat*kappa/sqrt(3) - 1i*sqrt( (sigmahat*kappa)^2/3 + 2i*kappa^2*lambdahat/3 );
alpha(2) = -1i*sigmahat*kappa/sqrt(3) + 1i*sqrt( (sigmahat*kappa)^2/3 + 2i*kappa^2*lambdahat/3 );
beta(1)  = -1i*sigmahat*kappa/sqrt(3) - 1i*sqrt( (sigmahat*kappa)^2/3 - 2i*kappa^2*lambdahat/3 );
beta(2)  = -1i*sigmahat*kappa/sqrt(3) + 1i*sqrt( (sigmahat*kappa)^2/3 - 2i*kappa^2*lambdahat/3 );

[~,imin] = min(real(alpha)); [~,imax] = max(real(alpha));
alphap = alpha(imin); alpham = alpha(imax);
[~,imin] = min(real(beta)); [~,imax] = max(real(beta));
betap = beta(imin); betam = beta(imax);

% first pass of guesses
Am = exp(alpham .* eta(1)); dAm = alpham*Am;
Bm = exp(betam .* eta(1)); dBm = betam*Bm;
Ap = exp(alphap*eta(end)); dAp = alphap*Ap;
Bp = exp(betap*eta(end)); dBp = betap*Bp;
AMP = 1;

    % marching forward
etam = eta(1:N/2); % left-half space
        % initialize
ym1 = zeros(length(etam),4);
ym2 = zeros(length(etam),4);
ym1(1,:) = AMP*[Am; dAm; 0; 0]; 
ym2(1,:) = AMP*[0; 0; Bm; dBm];
nm_orth = zeros(1,length(etam)); % index of normalization locations
Nm_orth = 0; % number of required orthonormalizations
for i = 2:length(etam)
    eta1 = etam(i-1); eta2 = eta(i-1)+1/2*deta; eta3 = eta2; eta4 = etam(i);
    
            % solution 1: B,dB = 0;
    y1 = ym1(i-1,:);
    k1 = deta*ddn(eta1,y1,lambdahat,sigmahat,kappa,deta).';
    y2 = ym1(i-1,:)+1/2*k1;    
    k2 = deta*ddn(eta2,y2,lambdahat,sigmahat,kappa,deta).';
    y3 = ym1(i-1,:)+1/2*k2;
    k3 = deta*ddn(eta3,y3,lambdahat,sigmahat,kappa,deta).';
    y4 = ym1(i-1,:) + k3;
    k4 = deta*ddn(eta4,y4,lambdahat,sigmahat,kappa,deta).';
    ym1(i,:) = ym1(i-1,:) + 1/6*( k1 + 2*k2 + 2*k3 + k4);            
%     [~,ym] = ode113(@ddn,[etam(i-1) etam(i)],ym1(i-1,:),[],l,s,Q0,deta);    
%     ym1(i,:) = ym(end,:);

            % solution 2: A,dA = 0;
    y1 = ym2(i-1,:);
    k1 = deta*ddn(eta1,y1,lambdahat,sigmahat,kappa,deta).';
    y2 = ym2(i-1,:)+1/2*k1;    
    k2 = deta*ddn(eta2,y2,lambdahat,sigmahat,kappa,deta).';
    y3 = ym2(i-1,:)+1/2*k2;
    k3 = deta*ddn(eta3,y3,lambdahat,sigmahat,kappa,deta).';
    y4 = ym2(i-1,:) + k3;
    k4 = deta*ddn(eta4,y4,lambdahat,sigmahat,kappa,deta).';
    ym2(i,:) = ym2(i-1,:) + 1/6*( k1 + 2*k2 + 2*k3 + k4);    
%     [~,ym] = ode113(@ddn,[etam(i-1) etam(i)],ym2(i-1,:),[],l,s,Q0,deta);
%     ym2(i,:) = ym(end,:);

        % criteria for orthonormalization
    YM = ym1(i,:) * ym2(i,:)'; % inner product takes c.c.
    YM1 = ym1(i,:) * ym1(i,:)';
    YM2 = ym2(i,:) * ym2(i,:)';
    a = acos( abs( YM / sqrt(YM1*YM2)) );
        % pass through orthonormalization process if necessary
    if a < a_max || i == length(etam);
            % turn on switch
        nm_orth(i) = i; Nm_orth = Nm_orth + 1;
        Y = [ym1(i,:).', ym2(i,:).'];
        [Z,P] = orthonorm(Y);
            % store matrices for reconstruction sweep
        eval(['Zm' num2str(Nm_orth) ' = Z;']);
        eval(['Pm' num2str(Nm_orth) ' = P;']);
            % take orthonormal vectors as initial conditions        
        ym1(i,:) = Z(:,1).';
        ym2(i,:) = Z(:,2).';
    end
end
    % marching backwords
etap = fliplr(eta(N/2:N)); % right-half space
        % initialize
yp1 = zeros(length(etap),4);
yp2 = zeros(length(etap),4);
yp1(1,:) = AMP*[Ap; dAp; 0; 0]; 
yp2(1,:) = AMP*[0; 0; Bp; dBp];
np_orth = zeros(1,length(etap)); % index of orthonormalization locations
Np_orth = 0; % number of required orthonormalizations
for i = 2:length(etap);
    eta1 = etap(i-1); eta2 = etap(i-1)-1/2*deta; eta3 = eta2; eta4 = etap(i);
    
            % solution 1: B,dB = 0;
    y1 = yp1(i-1,:);
    k1 = -deta*ddn(eta1,y1,lambdahat,sigmahat,kappa,deta).';
    y2 = yp1(i-1,:)+1/2*k1;    
    k2 = -deta*ddn(eta2,y2,lambdahat,sigmahat,kappa,deta).';
    y3 = yp1(i-1,:)+1/2*k2;
    k3 = -deta*ddn(eta3,y3,lambdahat,sigmahat,kappa,deta).';
    y4 = yp1(i-1,:) + k3;
    k4 = -deta*ddn(eta4,y4,lambdahat,sigmahat,kappa,deta).';
    yp1(i,:) = yp1(i-1,:) + 1/6*( k1 + 2*k2 + 2*k3 + k4);           
%     [~,yp] = ode113(@ddn,[etap(i-1) etap(i)],yp1(i-1,:),[],l,s,Q0,deta);
%     yp1(i,:) = yp(end,:);
            % solution 2: A,dA = 0;
    y1 = yp2(i-1,:);
    k1 = -deta*ddn(eta1,y1,lambdahat,sigmahat,kappa,deta).';
    y2 = yp2(i-1,:)+1/2*k1;    
    k2 = -deta*ddn(eta2,y2,lambdahat,sigmahat,kappa,deta).';
    y3 = yp2(i-1,:)+1/2*k2;
    k3 = -deta*ddn(eta3,y3,lambdahat,sigmahat,kappa,deta).';
    y4 = yp2(i-1,:) + k3;
    k4 = -deta*ddn(eta4,y4,lambdahat,sigmahat,kappa,deta).';
    yp2(i,:) = yp2(i-1,:) + 1/6*( k1 + 2*k2 + 2*k3 + k4);    
%     [~,yp] = ode113(@ddn,[etap(i-1) etap(i)],yp2(i-1,:),[],l,s,Q0,deta);
%     yp2(i,:) = yp(end,:);

        % criteria for orthonormalization
    YP = yp1(i,:) * yp2(i,:)';
    YP1 = yp1(i,:) * yp1(i,:)';
    YP2 = yp2(i,:) * yp2(i,:)';
    a = acos( abs( YP / sqrt(YP1*YP2)) );
        % pass through orthonormalization process if necessary
    if a < a_max || i == length(etap);
        np_orth(i) = i; Np_orth = Np_orth + 1;
        Y = [yp1(i,:).', yp2(i,:).'];
        [Z,P] = orthonorm(Y);
            % store matrices for reconstruction sweep
        eval(['Zp' num2str(Np_orth) ' = Z;']);
        eval(['Pp' num2str(Np_orth) ' = P;']);
            % take orthonormal vectors as initial conditions
        yp1(i,:) = Z(:,1).';
        yp2(i,:) = Z(:,2).';
    end
end



    % Reconstruction sweep
        % determine first set of constants for linear combination:
        % Solve matching conditions:
Ym1 = eval(['Zm' num2str(Nm_orth) '(:,1)']);
Ym2 = eval(['Zm' num2str(Nm_orth) '(:,2)']);
Yp1 = eval(['Zp' num2str(Np_orth) '(:,1)']);
Yp2 = eval(['Zp' num2str(Np_orth) '(:,2)']);

% Ym1 = ym1(end,:); Ym2 = ym2(end,:);
% Yp1 = yp2(end,:); Yp2 = yp2(end,:);

            % matching conditions
A = [Ym1(1), Ym2(1), -Yp1(1);
     Ym1(2), Ym2(2), -Yp1(2);
     Ym1(3), Ym2(3), -Yp1(3)];
B = [Yp2(1); Yp2(2); Yp2(3)];
C = A \ B;
C = [C;1];

        % error in matching condition
errM = ( C(1)*Ym1(4) + C(2)*Ym2(4) - C(3)*Yp1(4) - C(4)*Yp2(4)) / (C(4)*Yp2(4));

%%%%%%%%%%%%%%%%%%%%%%% Reconstruction process
        % Sweep backword from etam(end) to etam(1) for ym
            % set constants for linear combination
Cm = C(1:2);
            % get indices for reconstruction points
ym = zeros(length(etam),4);
ym(end,:) = eval(['Zm' num2str(Nm_orth)]) * Cm;
ind = find(nm_orth);
for i = 1:Nm_orth
        % refine constant for next set of orthornormalizations
    Cm = eval(['Pm' num2str(Nm_orth - (i-1))]) * Cm;
        % determine lower bound for current orthonormalization
    if i == Nm_orth;
        ii1 = 1;
    else ii1 = ind( Nm_orth-i );
    end
    ii2 = ind( Nm_orth-(i-1) ) - 1;
        % reconstruction process
    for ii = ii1 : ii2
        U = [ym1(ii,:).', ym2(ii,:).'];
        ym(ii,:) = U*Cm;
    end
end
        % Sweep forward from etap(end) to etap(1) for ym
            % set constants for linear combination
Cp = C(3:4);
            % get indices for reconstruction points
yp = zeros(length(etap),4);
yp(end,:) = eval(['Zp' num2str(Np_orth)]) * Cp;
ind = find(np_orth);
for i = 1:Np_orth
        % refine constant for next set of orthonormalizations
    Cp = eval(['Pp' num2str(Np_orth - (i-1))]) * Cp;
        % determine lower bound for current orthonormalization
    if i == Np_orth;
        ii1 = 1;
    else ii1 = ind( Np_orth-i );
    end
    ii2 = ind( Np_orth-(i-1) ) - 1;
        % reconstruction process
    for ii = ii1 : ii2
        U = [yp1(ii,:).', yp2(ii,:).'];
        yp(ii,:) = U * Cp;
    end
end  

%%%%%%%%%%%%% for debugging purposes, we plot
% figure(1)
% subplot(2,1,1)
% plot(etam,ym(:,1),etap,yp(:,1)); hold on
% plot(etam,ym(:,2),etap,yp(:,2));
% subplot(2,1,2)
% plot(etam,ym(:,3),etap,yp(:,3));
% plot(etam,ym(:,4),etap,yp(:,4));
% keyboard

end

