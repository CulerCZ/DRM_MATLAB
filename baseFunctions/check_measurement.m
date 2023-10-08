function [drp_measurement] = check_measurement(img_sample,drp_original,exp_para)
% this function is to check original DRPs acquired from measurement
% create date: Sep 22, 2021
% edit date: Sep 222, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------

figure('Name','demo_fig');
imshow(img_sample,'Border','tight');
colormap(jet)
[x,y] = ginput;
% press 'enter' to stop
nn = length(y);
y = fix(y);
x = fix(x);
close(findobj('type','figure','name','demo_fig'));

figure('Position',[200,200,200*(nn+1),200])
tiledlayout(1,nn+1,'TileSpacing','tight','Padding','compact')
ax = nexttile(1);
imshow(img_sample,'Border','tight')
colormap(ax,"jet")
hold on
scatter(x,y,108,'filled','w')
for ii = 1:nn
    text(x(ii)+5,y(ii)+5,int2str(ii),'FontSize',14)
end
hold off

drp_measurement = cell(nn,1);

for ii = 1:nn
    % DRP from measurement 
    nexttile(ii+1)
    x_pos = y(ii);
    y_pos = x(ii);
    drp_measurement{ii} = drp_original{x_pos,y_pos};
    DRPdisp(drp_measurement{ii},exp_para)
    
end
end

