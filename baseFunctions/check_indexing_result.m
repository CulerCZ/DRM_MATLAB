function [drp_measurement, drp_predicted, x, y] = check_indexing_result(EUmap,drp_original,exp_para,options)
% this function is to compare DRM indexing result between measured DRP and
% predicted DRP from indexing result.
% create date: Sep 6, 2021
% edit date: Sep 20, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------

arguments
    EUmap double
    drp_original cell
    exp_para struct
    options.plot_xtal (1,1) logical = false
end

euler = reshape(EUmap,size(EUmap,1)*size(EUmap,2),size(EUmap,3));
% get the color mapping of DRM measurement
cs = crystalSymmetry('cubic');
oM = ipfHSVKey(cs);
ori_drm_all = orientation.byEuler(euler(:,1)*degree,euler(:,2)*degree,euler(:,3)*degree,cs);
color_drm_all = oM.orientation2color(ori_drm_all);
color_drm_reg_all = reshape(color_drm_all,size(EUmap,1),size(EUmap,2),3);
clear euler oM ori_drm_all color_drm_all

figure('Name','demo_fig');
imshow(color_drm_reg_all,'Border','tight');
[x,y] = ginput;
% press 'enter' to stop
nn = length(y);
y = fix(y);
x = fix(x);
close(findobj('type','figure','name','demo_fig'));

figure('Position',[200,200,200*(nn+2),200*2.5])
tiledlayout(2,nn+2,'TileSpacing','tight','Padding','compact')
nexttile(1,[2,2])
imshow(color_drm_reg_all,'Border','tight')
hold on
scatter(x,y,72,'x','k')
for ii = 1:nn
    text(x(ii)+5,y(ii)+5,int2str(ii),'FontSize',14)
end
hold off

drp_measurement = cell(nn,1);
drp_predicted = cell(nn,1);

for ii = 1:nn
    % DRP from measurement 
    nexttile(ii+2)
    x_pos = y(ii);
    y_pos = x(ii);
    drp_measurement{ii} = drp_original{x_pos,y_pos};
    DRPdisp(DRP_norm(drp_measurement{ii}),exp_para)
    % DRP from prediction
    nexttile((ii+2+nn+2))
    eu_tmp = [EUmap(x_pos,y_pos,:)];
    if exist('drpTable_1')
        if fitQuality(x_pos,y_pos) > length(drpTable_1)
            drpsim_tmp = DRPsim_double(eu_tmp(1),eu_tmp(2),eu_tmp(3),fitting_para(1:4));
        else
            drpsim_tmp = DRPsim_double(eu_tmp(1),eu_tmp(2),eu_tmp(3),[4,0.01,16,4]);
        end
    else
        drpsim_tmp = DRPsim(eu_tmp(1),eu_tmp(2),eu_tmp(3),exp_para);
    end
    DRPdisp(drpsim_tmp,exp_para);
    drp_predicted{ii} = drpsim_tmp;
end

if options.plot_xtal
    plot crystal shape
    figure('Position',[200,200,200*nn,200])
    tiledlayout(1,nn,'TileSpacing','tight','Padding','compact')
    for ii = 1:nn
        % DRP from measurement 
        nexttile(ii)
        x_pos = y(ii);
        y_pos = x(ii);
        euler_angle = [EUmap(x_pos,y_pos,:)];
    %     cs = crystalSymmetry('cubic');
    %     ori_tmp = orientation.byEuler(eu_tmp(1)*degree, eu_tmp(2)*degree, eu_tmp(3)*degree,cs);
    %     cS = crystalShape.cube(cs);
    %     plot(ori_tmp * cS * 0.7);
    %     axis equal
    %     axis off
        eu1 = euler_angle(1);
        eu2 = euler_angle(2);
        eu3 = euler_angle(3);
        title(sprintf('%0.1f, %0.1f, %0.1f',[eu1, eu2, eu3]))

        drpsim_tmp = DRPsim(eu1,eu2,eu3,exp_para); % the matrix in form of th_num * ph_num
        C = [drpsim_tmp,drpsim_tmp(:,1)];
    %     theta = repmat([0:90/(exp_para.th_num-1):90]',1,exp_para.ph_num+1);
    %     phi = repmat([0:360/exp_para.ph_num:360],exp_para.th_num,1);
    %     x = cosd(theta).*cosd(phi);
    %     y = cosd(theta).*sind(phi);
    %     z = sind(theta);

    %     surf(x,y,z,C,'EdgeColor','none')
    %     axis equal
    %     set(gca,'visible','off')
    %     colormap('jet')
    %     alpha 0.45

        % Cube (for {100} faceting)
        r = .3;
        V1=[-1;  1; 1; -1; -1;  1; 1; -1;];
        V2=[-1; -1; 1;  1; -1; -1; 1;  1;];
        V3=[-1; -1;-1; -1;  1;  1; 1;  1;];
        F= [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8;];
        [THETA,PHI,R]=cart2sph(V1,V2,V3);
        R=r.*ones(size(V1(:,1)));
        [V1,V2,V3]=sph2cart(THETA,PHI,R);
        V=[V1 V2 V3];
        % apply the rotation
        roll = eu1*degree; pitch = eu2*degree; yaw = eu3*degree;
        dcm = angle2dcm(roll, pitch, yaw, 'ZXZ');
        V = V*dcm;
        patch('Faces',F,'Vertices',V,'FaceColor',[77 77 77]/255,'FaceAlpha',0.6,'EdgeColor','k',...
            'LineWidth',2); axis equal; grid on; hold on; view(3); 
        axis off
        view(0,20)
    end
end