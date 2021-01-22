function [K, M, shp] = alphaESM(xyz, alpha, mass, E)
% alphaESM: an Elastic Solid Model based on alpha shape 
% Author: Guang Song
% xyz: coordinates of atoms of the input structure
% alpha: alpha parameter for the alpha shape
% mass: mass of atoms of the input structure
% E: Young's modulus
% What to cite: 
% Guang Song, “Bridging Between Material Properties Of Proteins And The Underlying Molecular Interactions”, 2021
%
% Copyright (c) 2020 Guang Song. All Rights Reserved. 
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if nargin<4 E = 1.0; end
nu = 0.3; % default value for Poisson ratio.
% Now compute Lame parameters: lamdba, mu. 
lambda = nu*E/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
% compute alpha shape
shp = alphaShape(xyz, alpha);
elements = shp.alphaTriangulation;
% compute stiffness and mass matrices 
[K, volume]=stiffness_matrixP1_3D_elasticity( ...
    elements,xyz,lambda,mu);
K = (K+K')/2;
%K = full((K+K')/2); % make it fully symmetric
rho = sum(mass)/sum(volume);
M = mass_matrixP1_3D_elasticity(elements, volume*rho);
M = kron((M+M')/2, eye(3));
%M = full(M);
