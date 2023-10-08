function [drpDic, euDic, rotDic] = makeDRPdic(num,exp_para,options)
% the function is to generate DRP dictionary with corresponding labels
% Create date: Aug 17, 2021
% Edit date: Sep 20, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
arguments
    num {mustBeNumeric}
    exp_para struct
    options.engine string % choosing double or single
end

% th_max = exp_para.th_max;
% th_min = exp_para.th_min;
% th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

nn = 24*num;
randii = [rand(1,nn); rand(1,nn)];
thphii = [acosd(2*randii(1,:)-1)-90; 360*randii(2,:)];

vec = zeros(nn,3);
for ii = 1:nn
    tmp_thph = thphii(:,ii);
    vec(ii,:) = thph2vec(tmp_thph(1),tmp_thph(2));
end

quadr = and(and(vec(:,1)>=0,vec(:,2)>=0),vec(:,3)>=0);
ipf1 = and(vec(:,1)>=vec(:,2),vec(:,2)>=vec(:,3));
ipf2 = and(vec(:,2)>=vec(:,1),vec(:,1)>=vec(:,3));
good_idx = and(quadr,or(ipf1,ipf2));
goodvec = vec(good_idx==1,:);
goodnn = sum(good_idx);

drpDic = cell(goodnn,1);
euDic = zeros(goodnn,3);
rotDic = zeros(goodnn,1);

midpt = floor(ph_num/2);

tmp_vec_tmp = normr(rotate_facet(0,0,0,exp_para.faceting));
tmp_vec = tmp_vec_tmp(1,:);

parfor index = 1:goodnn
    vec1 = goodvec(index,:);
    vec2 = normr(cross(vec1,tmp_vec));
    vec3 = normr(cross(vec1,vec2));
    a_ij = [vec2', vec3', vec1'];
    
    euPh = acosd(a_ij(3,3));
    euph2 = atan2d(a_ij(1,3),a_ij(2,3));
    euph1 = atan2d(a_ij(3,1),a_ij(3,2));
    if options.engine == "single"
        tmpDRP = DRPsim(0,euPh,euph2,exp_para);
    else
        tmpDRP = DRPsim_double(0,euPh,euph2,exp_para);
    end
    colsum = sum(tmpDRP);
    [~,shift] = max(colsum);
    shift = 0 - shift;
    drp_uint8 = circshift(tmpDRP,shift,2);
    euDic(index,:) = [(360/ph_num)*shift,euPh,euph2];
    rotDic(index) = shift;
    % normalize the simulated DRP btw 0 and 1
    drpDic{index} = double(drp_uint8) / 255;
%     workbar(index/goodnn,'making DRP dictionary')
end
msgbox('DRP dictionary construction Completed!');