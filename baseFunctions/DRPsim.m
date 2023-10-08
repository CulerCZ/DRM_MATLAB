function [simDRP] = DRPsim(eu1,eu2,eu3,exp_para)
% the function is to simulate DRP based on Phong model.
% input Euler angle triplets are in Bunge convention.
% Edit date: Sep 14, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
arguments
    eu1 double
    eu2 double
    eu3 double
    exp_para struct
end

th_max = exp_para.th_max;
th_min = exp_para.th_min;
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;
faceting = exp_para.faceting;
i_Main = exp_para.fitting_para(1);
i_facet = exp_para.fitting_para(2);
sd_Main = exp_para.fitting_para(3);
sd_facet = exp_para.fitting_para(4);

rot_facet = normr(rotate_facet(eu1,eu2,eu3,faceting));
th_step = (th_max - th_min) / (th_num - 1);
ph_step = 360 / ph_num;

cauchy = @(p,x) p(1) ./ ((1+((x)./p(2)).^2));

th_range = th_min-th_step : th_step : ...
    (floor((90-th_min)/th_step)-1)*th_step+th_min;
th_DRP=repmat(transpose(th_range),1,ph_num);
ph_DRP=repmat(0:ph_step:360-ph_step,length(th_range),1);
simDRP=zeros(length(th_range),ph_num);
vec_DRP=zeros(3,length(th_range),ph_num);

% generate the vectors of incident/reflectance beams
for ii = 1:length(th_range)
    for jj = 1:ph_num
        tmp_vec = thph2vec(45+th_DRP(ii,jj)/2, ph_DRP(ii,jj));
        vec_DRP(:,ii,jj) = normr(tmp_vec);
    end
end

% major reflectance peak simulation
for ii = 1:length(rot_facet)
    ref_a1 = rot_facet(ii,:);
    ref = [0 0 -1] - 2 * dot([0 0 -1], ref_a1) * ref_a1;
    tmp_thph = vec2thph(ref);
    tmp_theta = tmp_thph(1);
    dPh=abs(ph_DRP - tmp_thph(2));
    dPh=abs(dPh - 360 * (dPh>180)); 
    peakDist=acosd(sind(th_DRP)*sind(tmp_theta)+cosd(th_DRP)*cosd(tmp_theta).*cosd(dPh));
    simDRP = max(simDRP, cauchy([i_Main, sd_Main], peakDist));
end

% great circle band simulation
for ii = 1:size(rot_facet,1)
    for jj = 1:size(rot_facet,1)
        if ii == jj
            continue
        end
        vec_1 = rot_facet(ii,:);
        vec_2 = rot_facet(jj,:);
        if all(vec_1 - vec_2 < 1e-3) || all(vec_1 + vec_2 < 1e-3)
            continue
        end
        gcnorm = normr(cross(vec_1, vec_2));
        peakDistb=zeros(size(simDRP));
        peakDista=zeros(size(simDRP));
        bandDist=zeros(size(simDRP));
        for mm = 1:length(th_range)
            for nn = 1:ph_num
                peakDista(mm,nn) = acosd(dot(vec_1, vec_DRP(:,mm,nn)));
                peakDistb(mm,nn) = acosd(dot(vec_2, vec_DRP(:,mm,nn)));
                bandDist(mm,nn) = asind(dot(gcnorm, vec_DRP(:,mm,nn)));
            end
        end
%         peakDist = bandDist; % + min(peakDista, peakDistb).^2.5;
        simDRP = max(simDRP, cauchy([i_facet, sd_facet], bandDist));
    end
end
    

simDRP=uint8(simDRP(2:th_num+1,:)*255);


end