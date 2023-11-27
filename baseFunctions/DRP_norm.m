function drp_norm = DRP_norm(drp)
% this funciton is to normalize a measured or simulated DRP to the scale of
% [0,1] in double type. The normlization will bring high-intensity features
% more substantial and weaken the minor features
% create date: Oct 18, 2021
% edit date: Oct 18, 2021
% Author: Chenyang ZHU @ NTU
arguments
    drp {mustBeNumeric}
%     exp_para struct
end

drp = double(drp);
upper = max(drp,[],'all');
lower = min(drp,[],'all');
gap = upper - lower;
drp_norm = uint8((drp - lower) / gap * 255);


end