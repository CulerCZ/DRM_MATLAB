function [all_rot] = rotate_facet(Aa,Bb,Cc,facet)
% to generate the equivalent orientation vectors with respect to (Aa,Bb,Cc)
% euler angle rotation direction pair.
% input x is a faceting vector x = [u v w], in order u >= v >= w.
% output all_rot is the whole vector combination.
% input Euler angle triplets are in Bunge convention.
% Edit date: Aug 27, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
if nargin < 4
    error('faceting information is not determined')
end

u = facet(1); v = facet(2); w = facet(3);

vec_eq = unique([u v w; -u v w; u -v w; u v -w;...
    u w v; -u w v; u -w v; u w -v; v w u; -v w u; v -w u; v w -u;...
    v u w; -v u w; v -u w; v u -w; w u v; -w u v; w -u v; w u -v;...
    w v u; -w v u; w -v u; w v -u],'row');
% this is for cubic system

N_vec = length(vec_eq);
% all_rot=EulerRotate(all11x,Aa,Bb,Cc);
all_rot = normr(EulerRotate(vec_eq,Aa,Bb,Cc));
for ii = 1:N_vec
    if all_rot(ii,3)<0
        all_rot(ii,:)=-all_rot(ii,:);
    end
end

end