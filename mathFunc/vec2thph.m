function [thph] = vec2thph(vec)
% Convert a vector [x,y,z] into theta and ph in degrees

normVec=normr(vec);
th=asind(normVec(3));
ph=atan2d(normVec(2),normVec(1))+360*(normVec(2)<0);

thph=[th ph];

end

