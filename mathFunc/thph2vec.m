function [vec] = thph2vec(th,ph)
% Convert th, ph in degrees into a vector [x y z]

nn = length(th);
vec=zeros(nn,3);
vec(:,1)=cosd(ph).*cosd(th);
vec(:,2)=sind(ph).*cosd(th);
vec(:,3)=sind(th);

end
