function [xyzr] = atoms2xyzr(xyz, atoms)
% compute a new xyzr matrix that contains the coordinates and radius of each atom 
% in the protein. 
% xyz: coordinates of atoms 
% atoms: element types of atoms
atomNames = 'HCNOS';
%radius2 = [1.2, 1.74, 1.54, 1.4, 1.8];
radius = [1.2, 2.0, 1.7, 1.6, 1.85]; % united atom radii
%mass = [1.0, 12, 14, 16, 32];

noH = find(atoms~='H');
xyz = xyz(noH,:);
atoms = atoms(noH,:);

% determine the radius for each atom
rz = zeros(size(atoms,1),1); 
for i=1:length(radius)
rz(find(atoms==atomNames(i)))=radius(i);
end
xyzr = [xyz, rz];