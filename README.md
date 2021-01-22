# sigmaESMrelease
The sigmaESM script and other MATLAB scripts for computing Young's moduli of proteins

The program assumes the following software, scripts, and data have been installed 
and added to the path using the addpath command of MATLAB.
(The program requires MATLAB 2020a or later.)
1. the MSMS software (available at http://mgltools.scripps.edu/downloads#msms)
2. the sbNMA release (available at https://github.com/htna/sbNMA-Matlab)
3. the sigmaESM release (available at https://github.com/gsongISU/sigmaESMrelease)
4. the dataset of the 18 proteins (list18.txt), their tertiary and quaternary structures, and the corresponding psfgen-generated psf/pdb files (available at https://github.com/gsongISU/sigmaESMrelease, under folder pdbDataset)
5. the two scripts for computing stiffness matrix 
	and mass matrix used in alphaESM (available at MATLAB file exchange https://www.mathworks.com/matlabcentral/fileexchange/27826-fast-fem-assembly-nodal-elements)

Run the MATLAB script sampleRun.m to compute the Young's moduli of the 18 proteins in the dataset and their interface regions (make sure to change the current folder to where the MSMS software is located.)


 
