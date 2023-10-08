function [drp_original] = igray2drp(igray_norm,phitheta,exp_para)
% the funtion is to convert igraynorm image stack into 4-d dataset in size
% of n1 * n2 * num_th * num_ph
% Edit date: Aug 27, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
[n1,n2,n3] = size(igray_norm);
if n3 ~= exp_para.th_num * exp_para.ph_num
    error('dimension is not match!')
end

drp_original = cell(n1,n2);
for ii = 1:n1
    for jj = 1:n2
        drp_original{ii,jj} = splitDRP(exp_para,igray_norm(ii,jj,:),phitheta);
    end
    workbar(ii/n1, sprintf('working on %d / %d rows',[ii n1]));
end

end