function [P, xyz] = rtbProjection(xyz, mass)
% the approach is to find the inertia. compute the principal axes. and then use them to determine directly translation and rotation. 
n = size(xyz, 1); % n: the number of atoms
if nargin == 1
  mass = ones(n,1);
end

M = sum(mass);
% find the mass center.
m3 = repmat(mass, 1, 3);
center = sum(xyz.*m3)/M;
xyz = xyz - center(ones(n, 1), :);

mwX = sqrt(m3).*xyz;
inertia = sum(sum(mwX.^2))*eye(3) - mwX'*mwX;
[V,D] = eig(inertia);
tV = V'; % tV: transpose of V. Columns of V are principal axes. 
for i=1:3
	trans{i} = tV(ones(n,1)*i, :); % the 3 translations are along principal axes 
end
P = zeros(n*3, 6);
for i=1:3
	rotate{i} = cross(trans{i}, xyz);
	temp = mat2vec(trans{i});
	P(:,i) = temp/norm(temp);
	temp = mat2vec(rotate{i});
	P(:,i+3) = temp/norm(temp);
end
m3 = mat2vec(sqrt(m3));
P = repmat(m3(:),1,size(P,2)).*P;
% now normalize columns of P
P = P*diag(1./normMat(P,1));

