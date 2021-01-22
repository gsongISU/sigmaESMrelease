% make sure to change the current folder to where the MSMS software is located)
charmm36prm; % load CHARMM36 parameters
temp = regexp(fileread('list18.txt'), '\r?\n', 'split');
list18 = vertcat(temp{:}); % read file list18.txt into a char array list18.

% the following block computes the young's modulus for each protein in the dataset
for i=1:18
  [cx1] = charmmPDBread([list18(i,:),'A_autopsf.pdb']);
  [atomTypes, bonds, angles, phis, imps] = psfread([list18(i,:),'A_autopsf.psf'], charmmAtomTypes);
  psfInfo = {bonds, angles, phis, imps, atomTypes};
  proteinE18(i) = proteinModulus0(cx1, psfInfo, charmmBonds, charmmAngles, charmmPhis, charmmImps, charmmNb, charmmMass);
end
% the following block computes the Young's modulus of each interface region of the 18 quaternary structures.
for i=1:18
  [cx1] = charmmPDBread([list18(i,:), '2_autopsf.pdb']);
  atomTypes = psfread([list18(i,:), '2_autopsf.psf'], charmmAtomTypes);
  chain1 = cx1(find(cx1(:,1)==1),3:5);
  chain2 = cx1(find(cx1(:,1)==2),3:5);
  [tempE, vol] = interfaceModulus0(chain1, chain2, atomTypes, charmmNb, charmmMass, cx1(:,6));
  interE(i)  = tempE;
  interVol(i)= vol;
end