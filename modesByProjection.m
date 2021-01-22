function [V, dia, bfr] = modesByProjection(Hess, mass, cx, nodeIdx)
% Hess: hessian matrix to be projected
% mass: masses of all the nodes
% cx: coordinates of the all the modes
% nodeIdx: a cell array that contains node indices in the project groups. The first two groups, representing the two protein regions, are treated as rigid 

% idx: indices of all atoms listed in the order of projection groups
idx = cell2mat(nodeIdx);
% idx3: the indices of atoms in modes: each atom is given 3 indices for its x, y, z
idx3 = [idx(:)*3-2, idx(:)*3-1, idx(:)*3]';
idx3 = idx3(:);
m3 = repmat(mass(:), 1, 3)';
m3sqrt = 1./sqrt(m3(:));
[i,j,s] = find(Hess);
N = size(Hess,1);
% now compute mass weighted hess
hessM = sparse(i,j,m3sqrt(i).*m3sqrt(j).*s,N,N);

% get projection matrix by groups
P = [];
for i=1:2
pp{i} = rtbProjectionMass(cx(nodeIdx{i},:), mass(nodeIdx{i}));
P = blkdiag(P, pp{i});
end
if length(nodeIdx) == 3
	P = blkdiag(P, eye(size(nodeIdx{3},1)*3));
end
% reorder P since idx is sorted by groups
reorderedP = zeros(size(P));
reorderedP(idx3,:) = P;

%P = blkdiag(pp{1},pp{2},pp{3});
Hp = reorderedP'*hessM*reorderedP;
Hp = (Hp+Hp')/2;
[Vp,D] = eig(full(Hp));
dia = diag(D);
% sort to make sure the eigenvalues are sorted in ascending order
[dia, idx] = sort(dia);
Vp = Vp(:, idx);
q = reorderedP*Vp(:,:);
V = repmat(m3sqrt, 1, size(q,2)).*q;
k = find(dia>1e-8, 1, 'first');
[bfr] = bfactorANM(V(:,k:end), dia(k:end));