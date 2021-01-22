function [hess, massMat, k_vdW, r, inter, intra] = sbNMA_inter(xyz1, xyz2, atomTypes, charmmNb, charmmMass) 
% sbNMA Hessian matrix for two interacting chains
% Author: Guang song
% xyz1: coordinates of chain1
% xyz2: coordiantes of chain2
% atomTypes: charmm atom types of both chain1 and chain2
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

% pairwise distance among all the atoms
n1 = size(xyz1,1);
n2 = size(xyz2,1);
n = n1+n2;
% vc: vdW contacts. 12: cutoff dist.
vc = sparse(n,n);
dist = pdist2(xyz1, xyz2);
dist = sparse(dist <= 12).*dist; % keep only contact distances <= 12
vc(1:n1,n1+1:n) = dist;
[vi, vj, r] = find(vc);
vc = [vi, vj];


[epsilon, r0] = computeKvdW(vc, [], atomTypes, charmmNb); 

tooClose = find(r<=r0);  
r(tooClose)=r0(tooClose); 
r6 = (r0./r).^6;
k_vdW = 12*epsilon.*(13*r6.^2-7*r6)./(r.^2);
k_vdW(find(k_vdW<0)) = 0; % zero all the negatives

cx_vdW = sparse(vc(:,1), vc(:,2), k_vdW, n, n);
cx12 =  (cx_vdW + cx_vdW');

hess=baseHessSparse([xyz1; xyz2], cx12);
id3 = mat2vec([atomTypes, atomTypes, atomTypes]);
massMat = charmmMass(id3+1); % mass matrix, +1: make it 1-based

function [epsilon, r0, containsH] = computeKvdW(vc, c14, atomTypes, charmmNb) % charmmAtomTypes)
charmmNb = charmmNb(find(charmmNb(:,1)>-1),:);
charmmNb = sortrows(charmmNb,1);
vc = atomTypes(vc) + 1; % make it 1-based
epsilon = sqrt(charmmNb(vc(:,1), 2).*charmmNb(vc(:,2), 2));
r0 = charmmNb(vc(:,1), 3) + charmmNb(vc(:,2), 3);

