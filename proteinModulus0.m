function [E, corr, massMat] = proteinModulus0(ubq, psfInfo, charmmBonds, charmmAngles, charmmPhis, charmmImps, charmmNb, charmmMass)
% this script computes the Young's modulus of a protein (intraE) using sbNMA and sigmaESM. 
% It calibrates the E value by requiring the magnitude of thermal vibrations 
% computed from sbNMA and sigmaESM be the same.
% ubq: a Nx7 matrix that contains information about each atom in the protein. col 1: chain number; col 2: atom index; cols 3-5: xyz; col 6: elment type (H, C, N, O, or S); col 7: backbone atoms ('A' for alpha carbon; 'C': C, 'N': N, 'O': O). 
% ubq is the direct output of charmmPDBreadread:
% e.g., [ubq] = charmmPDBread(['1p9gA_autopsf.pdb']);
% psfInfo: psf info for the structure. 
% charmm*: charmm parameters
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


% compute H and M matrices using sbNMA
[hess, massMat] = sbNMA_PSF(ubq(:,3:5), psfInfo{1}, psfInfo{2}, psfInfo{3}, psfInfo{4}, psfInfo{5}, charmmBonds, charmmAngles, charmmPhis, charmmImps, charmmNb, charmmMass);  
massSqrt = 1./sqrt(massMat(:));
hessM = massSqrt*massSqrt'.*hess;
%[Vm, Dm] = eig(hessM);
noH = find(ubq(:,6)~='H');
noH3 = find(mat2vec(repmat(ubq(:,6)~='H', 1, 3)));
invHess = pinv(full(hessM));
bfrHess = sum(vec2mat(diag(invHess)),2);
mass = massMat(1:3:end);
%bfrHess = bfrHess(noH); %commented out on 12/11/2020

xyzr = atoms2xyzr(ubq(:,3:5), ubq(:,6));
[verts,faces,volume] = msmsSurface0(xyzr);
E = 1.0;
[K] = sigmaESM(ubq(noH,3:5), verts, faces, volume, mass(noH), E);
massSqrtK = ones(size(K,1),1); % msms surface nodes given a mass of 1
massSqrtK(1:length(massSqrt(noH3))) = massSqrt(noH3); % the front portion of massSqrtK are atoms.
kM = massSqrtK*massSqrtK'.*K;
invK = pinv(full(kM));
bfrK = sum(vec2mat(diag(invK)),2);
massE = ones(size(K,1)/3, 1);
massE(1:length(noH))=mass(noH);
%bfrK = bfrK(1:length(noH));

try 
	corr = correlation(bfrHess, bfrK);
catch
	corr = 0;
end
bmagK = sum(bfrK)/sum(massE);
bmagHess = sum(bfrHess)/sum(mass);
E = bmagK/bmagHess*6.9;

