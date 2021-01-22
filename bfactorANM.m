function [bfr, weightedSF] = bfactorANM(V,dia,modes)
% Bfactor Calculation
% Author: Guang Song 
% Last update: 04/11/2017
% Last update: 02/2019
if nargin < 3
   modes = 1:size(V,2); % all the modes of V
end
Vsq = V(:,modes).^2;
% now compute squared fluctuations, or SF
SF = Vsq(1:3:end,:) + Vsq(2:3:end,:) + Vsq(3:3:end,:); % x^2 + y^2 + z^2
weightedSF = SF.*repmat((1./dia(modes))', size(SF,1), 1);
bfr = sum(weightedSF, 2);
%bfr = SF*(1./dia(modes)); % a N-by-1 vector
%bfr = sum(vec2mat(bfr), 2); % a N by 1 vector
