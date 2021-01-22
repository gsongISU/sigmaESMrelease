function [alpha] = findAlpha(cze, volume) 
% find an alpha value whose alphaShape volume is about the same as the given volume (<2% diff). 
% range: 2.0 to 5.2. intial guess: 3.6, the middle of the range
range = [2.0, 5.2];
while true 
	alpha = mean(range);
	shp = alphaShape(cze, alpha);
	elements = shp.alphaTriangulation;
	vol = sum(elementVolume(cze,elements));
	if abs(vol-volume)<0.02*volume % was 0.02
		break;
	elseif vol > volume
	  range(2) = mean(range);
	else
      range(1) = mean(range);
	end
end