function [volume] = elementVolume(nodes, elements)
% compute the volume of tetrahedral elements
% ref: https://www.mathworks.com/matlabcentral/answers/316266-how-i-can-calculate-area-and-volume-if-have-4-coordinates-like-x1-y1-x2-y2-x3-y3-x4-y4
%Volume = 1/6*abs(dot(cross(P2-P1,P3-P1),P4-P1));
for i=1:4
P{i} = nodes(elements(:,i),:);
end
volume = 1/6*abs(dot(cross(P{2}-P{1},P{3}-P{1}),P{4}-P{1}, 2));