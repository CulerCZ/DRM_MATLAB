function [eudata] = Euler2EBSD(EUmap, fname, non_index)
% convert the euler angle mapping into a readable file for mtex
% input EUmap is matrix of euler angle of every pixels in an image
% non_index_map: 0 for non-indexed; 1 for indexed.
% create date: Sep 06, 2021
% Edit date: Sep 14, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
arguments
    EUmap double
    fname string
    non_index logical
end

n1 = size(EUmap,1);
n2 = size(EUmap,2);
non_index_map = non_index;

eudata = zeros(n1*n2,6);
for ii = 1:n1
    for jj = 1:n2
        if non_index_map(ii,jj) == 0
            eudata((ii-1)*n2+jj,:) = ...
                [ii,jj,EUmap(ii,jj,1),EUmap(ii,jj,2),EUmap(ii,jj,3),0];
        else
            eudata((ii-1)*n2+jj,:) = ...
               [ii,jj,EUmap(ii,jj,1),EUmap(ii,jj,2),EUmap(ii,jj,3),1];
        end
    end
end

writematrix(eudata,fname,'Delimiter','tab')