function [E, interfaceVolume] = interfaceModulus0(chain1, chain2, atomTypes, charmmNb, charmmMass, atoms)
% chain1, chain2: the coordinates of chain1 and chain2
% atomTypes: the charmm atomTypes of both chain1 and chain2
% charmnNb: charmm non-bonded parameters
% charmmMass: charmm masses
% atoms: the standard atom types ('H', 'C', 'N', 'O', etc) of both chain1 and chain2
[hess, mass3] = sbNMA_inter2(chain1, chain2, atomTypes, charmmNb, charmmMass); % repmat: for 2 chains
mass = mass3(1:3:end);
n1 = size(chain1,1);
n2 = size(chain2,1);
whole = [chain1; chain2];
[V, dia, bfr] = modesByProjection(hess, mass, whole, {(1:n1)'; n1+(1:n2)'});
noH = find(atoms~='H'); 

xyzr = atoms2xyzr(whole, atoms);
[verts,faces,volume] = msmsSurface0(xyzr);
if volume == 0 % could not build msms surface
   E  = -1;
   interfaceVolume = 0;
   return;
end
E = 1.0;
[K, M, shp] = sigmaESM(whole(noH,:), verts, faces, volume, mass(noH), E);
elements = shp.alphaTriangulation;
noHatomsNumbers = [length(find(noH<=n1)), length(find(noH>n1))];
[nodes, interAll] = interfaceRegion(elements, noHatomsNumbers, shp.Points);
interfaceVolume = sum(elementVolume(shp.Points, elements(interAll,:)));

massSurfaceNodes = ones(size(shp.Points,1)-length(noH),1); % all MSMS surface nodes given a mass 1
massE = [mass(noH); massSurfaceNodes];
[V, dia, bfr2] = modesByProjection(K, massE, shp.Points, nodes);

% use mass weighted
E = (sum(bfr2.*massE)/sum(massE))/(sum(bfr.*mass)/sum(mass))*6.9


