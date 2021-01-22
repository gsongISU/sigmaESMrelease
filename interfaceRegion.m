function [nodes, interAll, intra1, intra2, inter, interface, surface] = interfaceRegion(elements, sz, cze)
% Determine the nodes and elements in the interface region. 
% elements: tetrahedral elements of the whole structure
% sz: the sizes (or chain lengths) of chain1 and chain2
% cze: the coordinates of atoms in chain1 and chain2, and surface nodes (in that order). 
% nodes: a cell array contains node indices of the two protein regions and the interface region. 
% intra1, intra2: elements of two protein chains, do not contain common nodes.
% inter: interface elements made of atoms from both protein chains
% interAll: indices of elements of the interface region
chainN = zeros(size(cze,1),1); % 2 protein chains + surface nodes
chainN(1:sz(1)) = 1;  % atoms in the first chain given a chainN of 1
chainN(sz(1)+1:sum(sz)) = -1; % second chain given a chainN of -1
chainId =  chainN(elements);
inter =  find((min(chainId,[], 2)== -1) & (max(chainId, [], 2)==1));
intra1 = find((max(chainId,[], 2)==1)   & min(chainId,[], 2)>= 0);
intra2 = find((min(chainId,[], 2)== -1) & max(chainId,[], 2)<= 0);
surface = find((min(chainId,[], 2)== 0) & (max(chainId, [], 2)==0));
e1 = elements(intra1,:);
e2 = elements(intra2,:);

common = intersect(e1(:), e2(:));
nodes1 = setdiff(unique(e1(:)), common);
nodes2 = setdiff(unique(e2(:)), common);
% find the rest that are not part of nodes1 or nodes2, which may be the same as common. 
rest = setdiff((1:size(cze,1))', [nodes1; nodes2]);
nodes = {nodes1; nodes2; rest};


intra = [intra1; intra2];
idx = [];
for i=1:length(common) 
   idx = [idx; find(sum(elements(intra,:)==common(i),2))];
end
idx = unique(idx);
interface = intra(idx);
intra1 = setdiff(intra1, interface);
intra2 = setdiff(intra2, interface);
%surface = setdiff(surface, interface2);

chainN(1:end) = 0;
chainN(nodes1) = 1;
chainN(nodes2) = -1;
chainId =  chainN(elements);
% interAll: include inter, interface, and part of the surface. 
interAll = find(abs(sum(chainId, 2))~=4);

