function DRPdisp(drp,exp_para,options)
% this function is to display a drp in polar coordinates
% drp is in size of th_num x ph_num
% Create date: Aug 27, 2021
% Edit date: Sep 14, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
arguments
    drp {mustBeNumeric}
    exp_para struct
    options.EdgeColor (1,1) string = 'none'
    options.colormap (1,1) string = 'jet'
    options.project (1,1) string = "stereo"
    options.cRange (1,2) double = [0 255]
    options.scaleBar (1,1) logical = 0
%     options.title (1,1) string = ""
end

if isfloat(drp(1))
    drp = uint8(drp.*255);
end

th_max = exp_para.th_max;
th_min = exp_para.th_min;
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

th_step = (th_max - th_min) / (th_num - 1);
ph_step = 360 / (ph_num - 1);

[x,y] = meshgrid(0:ph_step:360,th_min:th_step:th_max);
if options.project == "stereo"
    xx = cosd(y).*cosd(x)./(1+sind(y));
    yy = cosd(y).*sind(x)./(1+sind(y));
elseif options.project == "direct"
    xx = cosd(y).*cosd(x);
    yy = cosd(y).*sind(x);
end
h = pcolor(xx,yy,drp);
axis tight
axis equal
set(h, 'EdgeColor', options.EdgeColor);
ax = gca;
set(ax,'visible','off')
colormap(ax,options.colormap)
% title(h,options.title,'FontSize',16)
caxis(options.cRange)
if options.scaleBar
    colorbar('FontSize',16)
end

end