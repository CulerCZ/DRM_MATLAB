function [igray_norm] = bg_subtraction(igray_sample, igray_back, coeff)
% apply background subtraction
% Edit date: Aug 27, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
if nargin < 3
    coeff = 2;
end

n3 = size(igray_sample,3);
iback = zeros(size(igray_sample),'uint8');
for ii = 1:n3
    iback(:,:,ii) = wiener2(igray_back(:,:,ii),[7 7]);
    workbar(ii/n3,sprintf('removing noise %d / %d',[ii,n3]));
end
igray_norm = zeros(size(igray_sample),'uint8');
for ii = 1:n3
    igray_norm(:,:,ii) = uint8(double(igray_sample(:,:,ii))./...
        double(iback(:,:,ii)) / coeff * 255);
    
    workbar(ii/n3,sprintf('generating igray_norm %d / %d',[ii n3]));
end

igray_norm(igray_back<10) = 0;

end