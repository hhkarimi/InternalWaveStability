function [ Q ] = beam_profile( x )
%BEAM_PROFILE beam profile of peak amplitude 1
%   input:  x = column vector of spatial coordinate
%   output: Q = column vector complex beam profile at x


%     % standing vs progressive beam
% k = linspace(0,20,1000); % progressive beam
% % k = linspace(-10,10, 2000); % standing beam
% 
% f = 0.1;
% 
%     % Taking Gaussian profile in wavenumber space
%         % proper normalization
% b = 1/8;
% a = b * sqrt(2/pi) * exp(1/2);
% Fk = a * exp( -k.^2 * b);
% 
%     % integrate over wavenumbers for each x
% Q = zeros(length(x),1);
% for i = 1:length(x);
%     I = Fk .* exp(1i * k * x(i));
%     Q(i) = trapz(k,I);
% end

% or simply Gaussian
Q = 1/2 * exp(-x.^2);


end