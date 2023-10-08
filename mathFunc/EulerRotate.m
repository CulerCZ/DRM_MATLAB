function [U2]= EulerRotate(U1,Aa,Bb,Cc)
% takes a row vector and 3 euler angles in degree
% returns the rotation of that row vector through the three angles
% Edit date: Aug 28, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
Aa=double(Aa);
Bb=double(Bb);
Cc=double(Cc);

Rxb=[1 0 0; 0 cosd(Bb) sind(Bb); 0 -sind(Bb) cosd(Bb)];
Rza=[cosd(Aa) sind(Aa) 0; -sind(Aa) cosd(Aa) 0; 0 0 1];
Rzc=[cosd(Cc) sind(Cc) 0; -sind(Cc) cosd(Cc) 0; 0 0 1];

OoO=Rzc*Rxb*Rza;

U2=U1*OoO;

% coordinate transformation from crystal reference frame into sample
% reference frame.

end