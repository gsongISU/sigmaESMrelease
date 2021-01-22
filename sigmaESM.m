function [K, M, shp] = sigmaESM(xyz, verts, faces, volume, mass, E)
% simgaESM: ESM based on molecular surface 
% Author: Guang Song
% K, M: stiffness and mass matrices of the elastic solid model
% shp: elastic solid model represented by an alpha shape object 
% xyz: atomic coordinates of the input structure
% verts, faces, volume: vertices and faces of MSMS surface and its enclosed volume 
% mass: masses of atoms of the input structure
% E: Young's modulus of the solid model. Set to 1 (unit: kcal/mol/A^3) by dafault.
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

if nargin<6 E = 1.0; end
init_alpha = 3.6; % 2.1 atom radius + 1.5 (probe radius)
shp0 = alphaShape(xyz, init_alpha); 
[pface] = boundaryFacets(shp0);
m = size(pface,1); % m: # of faces on the initial alpha shape
% simplify the molecular surface to have m faces  
[rF, rV] = reducepatch(faces, verts, m); 
xyz = [xyz; rV]; % append reduced surface vertices rV to xyz
% tune the alpha so that the alpha shape has the given volume
alpha = findAlpha(xyz, volume);
[K, M, shp] = alphaESM(xyz, alpha, mass, E);