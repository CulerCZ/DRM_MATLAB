%% DRM measurement indexing engine
exp_para.th_max = 65;
exp_para.th_min = 8;
exp_para.th_num = 20;
exp_para.ph_num = 72;
exp_para.ph_min = 0;
exp_para.ph_max = 355;
exp_para.faceting = [1 1 1];
% exp_para.fitting_para = [1, 0.6, 15, 8, 0.8, 8];
% for nickel, the fitting parameters:
exp_para.fitting_para = [1 0.5 15 6 0.8 8];
%% load sample and background dataset
scaleCoeff = 0.5;
[igray_sample_top, phitheta, pos_top, img_sample_top] = drp_loader( ...
    exp_para,zeros(1,4),format='jpg',scale=scaleCoeff);
[igray_sample_bot, ~, pos_bot, img_sample_bot] = drp_loader( ...
    exp_para,zeros(1,4),format='jpg',scale=scaleCoeff);
%%
[movingPoints, refPoints] = cpselect(img_sample_bot,img_sample_top,'Wait',true);
tform_2surfaces = fitgeotrans(movingPoints,refPoints,'affine');
output_region = imref2d(size(img_sample_top));
img_back_trans = imwarp(img_sample_bot,tform_2surfaces,'nearest','OutputView',output_region);
figure(11), imshowpair(img_sample_top,img_back_trans,'montage')
% exportgraphics(gcf,fullfile(saveFolder,"region_registering_r.tiff"),Resolution=300)
close(11)
igray_sample_bot_trans = igray_sample_top;
for ii = 1:size(igray_sample_bot,3)
    back_temp = igray_sample_bot(:,:,ii);
    igray_sample_bot_trans(:,:,ii) = imwarp(back_temp,tform_2surfaces,'nearest','OutputView',output_region);
    if mod(ii,100)==0
        fprintf("finish %d / %d ...\n",[ii,size(igray_sample_bot,3)]);
    end
end
clear back_temp
%% select gauge range and crop the igray dataset
figure(101)
imshow(igray_sample_top(:,:,end) * 2,'Border','tight')
roi = drawrectangle;
pos_crop = roi.Position;
close 101
temp_img_crop = imcrop(igray_sample_top(:,:,end),pos_crop);
[n1,n2] = size(temp_img_crop);
igray_f = zeros(n1,n2,size(igray_sample_top,3),'uint8');
igray_b = igray_f;
for ii = 1:size(igray_sample_top,3)
    igray_f(:,:,ii) = imcrop(igray_sample_top(:,:,ii),pos_crop);
    igray_b(:,:,ii) = imcrop(igray_sample_bot_trans(:,:,ii),pos_crop);
    if mod(ii,100)==0
        fprintf("finish %d / %d ...\n",[ii,size(igray_sample_bot,3)]);
    end
end

%% background subtratction applied
% convert into drps' stack
% drp_original is in size of [n1 x n2 x th_num x ph_num]
drp_original_top = igray2drp(igray_f,phitheta,exp_para);
drp_original_bot = igray2drp(igray_b,phitheta,exp_para);
% uisave('igray_norm','igray_norm')

%% quick check DRP and show DRP
drp_measurement = check_measurement(img_sample_bot,drp_original_bot,exp_para); 
%% show sim DRP
% ori_temp = grain_left.meanOrientation;
% eu1 = ori_temp.phi1/degree;
% eu2 = ori_temp.Phi/degree;
% eu3 = ori_temp.phi2/degree; bv
figure('Position',[200,200,200,200])
DRPdisp(DRPsim(49.5,45,30,exp_para),exp_para);

%% generate DRP dictionary
num_dic = 10000;
[drpDic, euDic, rotDic] = makeDRPdic(num_dic,exp_para,engine="single");
% [drpDic_2, euDic_2, rotDic_2] = makeDRPdic(num_dic,exp_para,engine="double");
%
clear num_dic

%% training an autoencoder
hiddenSize1 = 100;
AE_DRM = trainAutoencoder(drpDic,hiddenSize1, ...
    'MaxEpochs',200, ...
    'L2WeightRegularization',0.001, ...
    'SparsityRegularization',4, ...
    'SparsityProportion',0.10, ...
    'ScaleData', false, ...
    'UseGPU', false);
clear hiddenSize1
%%
% start indexing engine and obtain final result
index_result_top = IndexingEngine(drp_original_top,AE_DRM,exp_para,drpDic,euDic,rotDic);
index_result_bot = IndexingEngine(drp_original_bot,AE_DRM,exp_para,drpDic,euDic,rotDic);
figure, imshow(plot_ipf_map(index_result_top.EUmap),'Border','tight')
% exportgraphics(gcf,"D:\DRM\rawData\groupB_Sample04_top\ipf_top.tif",Resolution=300)
% close all
figure, imshow(plot_ipf_map(index_result_bot.EUmap),'Border','tight')
% exportgraphics(gcf,"D:\DRM\rawData\groupB_Sample04_top\ipf_bot.tif",Resolution=300)
% close all
%% check indexing results
[drp_measurement, drp_predicted, x, y] = check_indexing_result(index_result.EUmap,drp_original,exp_para);
%%
drp_out = drp_encode(AE_DRM,drp_measurement,exp_para);
drp_tmp = drp_encode(AE_DRM,drp_predicted,exp_para);

%% check indexing outcome
check_indexing_result(indexResult_top.euMap,drp_original_top,exp_para);


%% create boundary information
[boundary_front,boundary_map_front,grain_map_front,icanny_front] = calc_gb(drp_f,exp_para);
[boundary_back,boundary_map_back,grain_map_back,icanny_back] = calc_gb(drp_b,exp_para);
figure, imshow(grain_map_front)
figure, imshow(grain_map_back)

%% select boundaries
figure, imshowpair(boundary_map_front>0,boundary_map_back>0);
hold on
% imshow(boundary_map_back);
gb_corr_selected = gb_corr(gb_corr(:,1) > 0 & gb_corr(:,2) > 0,:);
for ii = 1:length(gb_corr_selected)
    idx_1 = gb_corr_selected(ii,1);
    idx_2 = gb_corr_selected(ii,2);
%     if idx_2 == 0
%         continue
%     end
    line([centroid_f(idx_1,2),centroid_b(idx_2,2)],...
        [centroid_f(idx_1,1),centroid_b(idx_2,1)],'Color','w','LineWidth',3)
end
hold off
% exportgraphics(gcf,fullfile(saveFolder,"GBpaired.tif"),Resolution=300)
[GBC,EUmap_front,EUmap_back] = calc_gbc(gb_corr_selected,boundary_front,boundary_back,index_result_f,index_result_b);

%%
figure, hold on
nn = length(GBC);
cmap = colormap(hot);
for ii = 1:length(GBC)
    len_front = length(GBC{ii,1});
    len_back = length(GBC{ii,2});
    len_curve = max(len_front,len_back);
    % front part
%     [xx, yy] = interp_GB_trace(GBC{ii,1}(:,1), GBC{ii,1}(:,2), len_curve, smooth=false);
%     GBC{ii,6} = [xx,yy];
% %     line(xx,yy,'Color','r','LineWidth',2)
% %     scatter(xx,yy,2,'filled','MarkerFaceColor',ii/nn*[1 0 0]+(1-ii/nn)*[0 0 1])
%     cfrac = min(1,GBC{ii,5} / 62);
%     scatter(yy,xx,3,'filled','MarkerFaceColor',cmap(fix(cfrac*length(cmap)),:))
    [xx, yy] = interp_GB_trace(GBC{ii,2}(:,1), GBC{ii,2}(:,2), len_curve, smooth=false);
    GBC{ii,7} = [xx,yy];
%     line(xx,yy,'Color','b','LineWidth',2)
    cfrac = min(1,GBC{ii,5} / 62);
    scatter(yy,xx,3,'filled','MarkerFaceColor',cmap(fix(cfrac*length(cmap)),:))
% %     workbar(ii/length(GBC),sprintf('processing %d / %d GB...',[ii length(GBC)]));
end
axis equal
scale_ratio = 4; % mm in total
xmax = size(boundary_front.grainIdxMap,1);
ylim([0 xmax])
ymax = size(boundary_front.grainIdxMap,2);
xlim([0 ymax])
box on
set(gca,'LineWidth',1.5)
% colormap(hot)
% colorbar('Ticks',[0 1],'TickLabels',{'0','62'})
% x_len = 41.5;
% map_size = size(boundary_front.grainIdxMap);
% [x_ticks, y_ticks, x_ticklabels, y_ticklabels] = generate_tick(x_len, map_size, interval=10);
xticks([])
% xticklabels(split(x_ticklabels))
yticks([])


%% plot grain boundary traces in 3-d
thickpxlratio = 15;
grainNum = 3;
pixelBot = GBC{grainNum,6};
[~,seq] = sort(pixelBot(:,1),'ascend');
pixelBot = pixelBot(seq,:);
pixelTop = GBC{grainNum,7};
[~,seq] = sort(pixelTop(:,1),'ascend');
pixelTop = pixelTop(seq,:);
pixelTilt = GBC{grainNum,8};
pixelNum = length(pixelBot);
% for ii = 1:pixelNum-1
% 
% end
xrange = range([pixelTop(:,1);pixelBot(:,1)]);
yrange = range([pixelTop(:,2);pixelTop(:,2)]);
xmin = min([pixelTop(:,1);pixelBot(:,1)]);
ymin = min([pixelTop(:,2);pixelTop(:,2)]);
voxels = zeros(xrange+1,yrange+1,thickpxlratio+1);

% define a voxel on grain boundary plane or not
for iix = 1:xrange+1
    for iiy = 1:yrange+1
        for iiz = 1:size(voxels,3)
            ix = iix-1;
            iy = iiy-1;
            iz = iiz-1;
            voxelCoord = [xmin+ix, ymin+iy, 0+iz];
            disList = zeros(1,pixelNum);
            for ii = 1:pixelNum
                lineDir = [pixelTop(ii,:),thickpxlratio] - [pixelBot(ii,:),0];
                toPointVec = voxelCoord - [pixelBot(ii,:),0];
                proj_temp = (lineDir*toPointVec')/norm(lineDir)^2 * lineDir;
                disList(ii) = norm(toPointVec - proj_temp);
            end
            voxels(iix,iiy,iiz) = min(disList);
        end
    end
end

% illustration of the boundary traces
% figure,
% tiledlayout(4,4,'TileSpacing','tight','Padding','tight')
% for ii = 1:16
%     nexttile(ii)
%     imshow(voxels(:,:,ii)<=1)
% end
% exportgraphics(gcf,sprintf("/Users/chenyangzhu/Desktop/AlGBDRM/GBFigs/GBsection_%d.tif",grainNum),'Resolution',300)
%
[gbx, gby, gbz] = ind2sub(size(voxels),find(voxels<=1));
gbVoxNum = length(gbx);
figure, hold on
for ii = 1:gbVoxNum
    plotCube([gbx(ii),gby(ii),gbz(ii)],1,color=[173, 216, 230]/256,facealpha=0.717)
end
axis equal
xlim([-1 xrange+1]) 
zlim([-1 thickpxlratio+1])
ylim([-1 yrange+1])
set(gca,'visible','off')

plotBox([0 0 0],size(voxels,1),size(voxels,2),size(voxels,3))
view(-110,30)
% exportgraphics(gcf,sprintf("/Users/chenyangzhu/Desktop/AlGBDRM/GBFigs/GB3d_%d.tif",grainNum),'Resolution',300)
exportgraphics(gcf,fullfile(saveFolder,sprintf("GB_sample_3d_%d.tif",grainNum)),Resolution=300)
%
figure
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile(1)
plot_boundary(boundary_front, gb_corr_selected(grainNum,1))
nexttile(2)
plot_boundary(boundary_back, gb_corr_selected(grainNum,2))
% exportgraphics(gcf,sprintf("/Users/chenyangzhu/Desktop/AlGBDRM/GBFigs/GB2d_%d.tif",grainNum),'Resolution',300)
exportgraphics(gcf,fullfile(saveFolder,sprintf("GB_sample_2d_%d.tif",grainNum)),Resolution=300)


%% functions used in tensile_sample_code
function drpLib = DRPLibGenerator(ang_res, exp_para, options)
% this function requires MTEX package for generation of orientation library
    arguments
        ang_res (1,1) double 
        exp_para (1,1) struct
        options.verbose (1,1) logical = 1
        options.crystalSymmetry (1,:) string = 'cubic'
    end
    cs = crystalSymmetry(options.crystalSymmetry);
    ori = equispacedSO3Grid(cs,'resolution',ang_res);
    nn = length(ori.phi1);
    drpLib.drpDic = cell(nn,1);
    drpLib.eulerDic = zeros(nn,3);
    for ii = 1:nn
        eu = [ori(ii).phi1, ori(ii).Phi, ori(ii).phi2]./degree;
        drpLib.drpDic{ii} = DRPsim(eu(1),eu(2),eu(3),exp_para);
        drpLib.eulerDic(ii,:) = eu;  % in degree
        if options.verbose
            workbar(ii/nn,sprintf("processing %d / %d DRPs",[ii nn]));
        end
    end
end


function indexResult = DirectDIEngine(drpM, drpLib, options)
    arguments
        drpM
        drpLib
        options.K (1,1) double = 1
    end
    [n1,n2] = size(drpM);
    EUmap = zeros(3,n1,n2);
    drplist_s = zeros(length(drpLib.drpDic), numel(drpM{1,1}));
    for ii = 1:length(drpLib.drpDic)
        drplist_s(ii,:) = double(reshape(drpLib.drpDic{ii},1,[]))/256;
    end
    for ii = 1: n1
        drplist_m = zeros(n2,numel(drpM{1,1}));
        for jj = 1:n2
            drplist_m(jj,:) = double(reshape(drpM{ii,jj},1,[]))/256;
        end

        [Idx, D] = knnsearch(drplist_s, drplist_m, K=options.K);

        EUmap(:,ii,:) = drpLib.eulerDic(Idx(:,1),:)';
        indexResult.idxMap(ii,:) = Idx(:,1);
        workbar(ii/n1, sprintf("processing %d / %d lines...",[ii n1]));
    end
    indexResult.euMap = permute(EUmap, [2,3,1]);
end


function [boundary,boundary_map,grain_map,icanny] = calc_gb(drp_original,exp_para,options)
    arguments
        drp_original cell
        exp_para struct
        options.conn (1,1) double = 4
        options.canny_filter (1,1) double = 10
    end
    % first step is to calculate grain map from drp_original
    th_num = exp_para.th_num;
    ph_num = exp_para.ph_num;
    [n1,n2] = size(drp_original);
    thr_canny = [0.2 0.4];
    sigma = 4;
    icanny = zeros(n1,n2);
    for ii = 1:th_num
        for jj = 1:ph_num
            layer_fig = extract_fig(drp_original,exp_para,[ii,jj]);
            img_gb = edge(layer_fig,'canny',thr_canny,sigma);
            icanny = icanny + img_gb;
        end
        workbar(ii/th_num, sprintf('processing %d / %d ...',[ii th_num]));
    end
    boundary_map = icanny > options.canny_filter;
    grain_skel = bwskel(boundary_map);
    grain_map = set_boundary(~grain_skel,width=5);
    cc = bwconncomp(grain_map,4);
    grain_index_map = bwlabel(grain_map,4); % grain label from 0 to grain_num
    grain_info = cell(cc.NumObjects,1);
    % iterate over all filtered grains
    for ii = 1:cc.NumObjects
        % use peripheral pixels of the grains, the intersect of which are
        % grain boundary between two adjacent grains
        mask = [1,1,1;1,8,1;1,1,1];
        gg = filter2(mask,grain_index_map == ii,'same');
        bb = gg>0 & gg<8; % each grain with its boundary
%         tmp_boundary_pos = bwboundaries(bb,8); % use 8-connect
        cc_tmp = bwconncomp(bb,4);
        % store the label of pixels on grain boundary
        grain_info{ii} = cc_tmp.PixelIdxList{1};
        workbar(ii/cc.NumObjects,sprintf('processing grains %d / %d',[ii cc.NumObjects]))
    end
    
    boundary = struct(...
        'NumOfBoundaries', 0, ...
        'PixelIdxList', [], ...
        'PixelList', [], ...
        'NumPixels', [], ...
        'grainAdjacent', [], ...
        'tripleJunction', [], ...
        'grainIdxMap', grain_index_map);
    idx = 1;
    boundary_pixels = [];
    for ii = 1:cc.NumObjects
        for jj = ii+1:cc.NumObjects
            p1 = grain_info{ii}';
            p2 = grain_info{jj}';
            if length(p1) < 10 || length(p2) < 10
                continue
            else
                grain_pixel_idx = intersect(p1,p2);
                if ~isempty(grain_pixel_idx)
                    if ~isempty(intersect(grain_pixel_idx,boundary_pixels))
                        boundary.tripleJunction = ...
                            [boundary.tripleJunction, intersect(grain_pixel_idx,boundary_pixels)];
        %                 grain_pixel_idx(grain_pixel_idx == intersect(grain_pixel_idx,boundary_pixels)) = [];
                    end
                    boundary_pixels = unique([boundary_pixels, grain_pixel_idx]);
                    boundary.NumPixels(idx) = numel(grain_pixel_idx);
                    boundary.PixelIdxList{idx} = grain_pixel_idx;
                    boundary.grainAdjacent{idx} = [ii,jj];
                    idx = idx + 1;
                end
            end
        end
        workbar(ii/cc.NumObjects,sprintf('processing grains %d / %d',[ii cc.NumObjects]))
    end
    boundary.NumOfBoundaries = idx - 1;
    boundary_map = zeros(n1,n2);
    num_of_boundary = length(boundary.PixelIdxList);
    for ii = 1:num_of_boundary
        boundary_map(boundary.PixelIdxList{ii}) = ii/(num_of_boundary+1);
        [r,c] = ind2sub(size(drp_original),boundary.PixelIdxList{ii});
        boundary.PixelList{ii} = zeros(length(r),2);
        boundary.PixelList{ii}(:,1) = r';
        boundary.PixelList{ii}(:,2) = c';
    end
end

function layer_fig = extract_fig(drp_original,exp_para,pos,options)
    % the function is to extract layered figures from drp_original cell
    % dataset. The input position is the (phi,theta) combination so the
    % extracted figure is from observation angle (phi,theta).
    arguments
        drp_original cell
        exp_para struct
        pos (1,2) double = [1 1]
        options.dtype (1,:) string = 'jpg'
    end
    th_num = exp_para.th_num;
    ph_num = exp_para.ph_num;
    [n1,n2] = size(drp_original);
    layer_fig = zeros(n1,n2,'uint8');
    if pos(1) <= 0 || pos(1) > th_num || pos(2) <= 0 || pos(2) > ph_num
        error('position out of range!')
    end
    for ii = 1:n1
        for jj = 1:n2
            layer_fig(ii,jj) = drp_original{ii,jj}(pos(1),pos(2));
        end
    end
end

function grain_map_r = set_boundary(grain_map,options)
% set grain map boundaries to be zero, which are not grains anymore`
    arguments
        grain_map logical
        options.width (1,1) double = 5
    end
    grain_map_r = grain_map;
    [n1,n2] = size(grain_map);
    w = options.width;
    grain_map_r(:,1:w) = 0;
    grain_map_r(1:w,:) = 0;
    grain_map_r(:,n2-5:n2) = 0;
    grain_map_r(n1-5:n1,:) = 0;
end

function [gb_corr, centroid_1, centroid_2] = pair_gbn(boundary_1, boundary_2, options)
% this function is to pair two sides of grain boundary network
% automatically. The output gb_corr will be the correlation matrix to
% connect grain boundaries that are successfully paired.
    arguments
        boundary_1 (1,1) struct
        boundary_2 (1,1) struct
        options.length_thres (1,1) double = 10 % unit in pixels
        options.n_out (1,1) int8 = 1
    end
    n1 = boundary_1.NumOfBoundaries;
    n2 = boundary_2.NumOfBoundaries;
    centroid_1 = zeros(n1,3);
    centroid_2 = zeros(n2,3);
    for ii = 1:n1
        centroid_1(ii,1:2) = mean(boundary_1.PixelList{ii});
        centroid_1(ii,3) = find_slope(centroid_1(ii,1:2));
    end
    for ii = 1:n2
        centroid_2(ii,1:2) = mean(boundary_2.PixelList{ii});
        centroid_2(ii,3) = find_slope(centroid_2(ii,1:2));
    end
    gb_corr = zeros(n1,options.n_out+1); % store correlation matrices...not sure yet
    for ii = 1:n1
        x_1 = centroid_1(ii,1);
        y_1 = centroid_1(ii,2);
        len_1 = boundary_1.NumPixels(ii);
%         len_2 = boundary_2.NumPixels(ii);
        if len_1 < options.length_thres
            continue
        else
            rel_pos = [x_1 y_1] - centroid_2(:,1:2);
            dis_tmp = sum(rel_pos.*rel_pos, 2);
            [~, idx] = mink(dis_tmp,5);
%             idx(abs(atand(centroid_2(idx,3))-atand(centroid_1(ii,3)))>20) = 0;
            idx(boundary_2.NumPixels(idx) < options.length_thres) = 0;
            gb_corr(ii,:) = [ii idx(1:options.n_out)'];
        end
    end
    
end

function slope = find_slope(pos_cord)
    arguments
        pos_cord (:,2) double
    end
    slope = pos_cord(:,1) \ pos_cord(:,2);
end

function [GBC,EUmap_front,EUmap_back] = calc_gbc(gb_corr,boundary_front,boundary_back,index_result_front,index_result_back,options)
% the functino is to calculate grain boundary characteristics with
% correlation matrix and indexing results
    arguments
        gb_corr double
        boundary_front struct
        boundary_back struct
        index_result_front struct
        index_result_back struct
        options.smooth (1,1) logical = true
        options.crystalsymmetry string = 'cubic'
%         options.crop (1,1) logical = false
        options.range (1,4) double
    end
    % the columns of output GBC will be trace_front | trace_back | ori_1 |
    % ori_2 | disorientation angle (for now, need to be updated)
    
    num_gb = length(gb_corr);
    GBC = cell(num_gb,4);
    try
        x1 = options.range(1);
        x2 = options.range(2);
        y1 = options.range(3);
        y2 = options.range(4);
        EUmap_tmp_front = permute(index_result_front.EUmap(x1:x2,y1:y2,:),[3,1,2]);
        EUmap_tmp_back = permute(index_result_back.EUmap(x1:x2,y1:y2,:),[3,1,2]);
    catch
        EUmap_tmp_front = permute(index_result_front.EUmap,[3,1,2]);
        EUmap_tmp_back = permute(index_result_back.EUmap,[3,1,2]);
    end
    % preassign average grain orientation map
    EUmap_front = zeros(size(EUmap_tmp_front));
    EUmap_back = zeros(size(EUmap_front));
    
    for idx = 1:num_gb
        idx_front = gb_corr(idx,1);
        idx_back = gb_corr(idx,2);
        GBC{idx,1} = boundary_front.PixelList{idx_front};
        GBC{idx,2} = boundary_back.PixelList{idx_back};
        % grain information
        grain_1_front = boundary_front.grainIdxMap==boundary_front.grainAdjacent{idx_front}(1);
        euler_1_front = EUmap_tmp_front(:,grain_1_front);
        grain_1_back = boundary_back.grainIdxMap==boundary_back.grainAdjacent{idx_back}(1);
        euler_1_back = EUmap_tmp_back(:,grain_1_back);
        grain_2_front = boundary_front.grainIdxMap==boundary_front.grainAdjacent{idx_front}(2);
        euler_2_front = EUmap_tmp_front(:,grain_2_front);
        grain_2_back = boundary_back.grainIdxMap==boundary_back.grainAdjacent{idx_back}(2);
        euler_2_back = EUmap_tmp_back(:,grain_2_back);
        if options.smooth
            cs = crystalSymmetry(options.crystalsymmetry);
            ori_1_front = orientation.byEuler(euler_1_front'.*degree,cs);
            ori_1_back = orientation.byEuler(euler_1_back'.*degree,cs);
            ori_2_front = orientation.byEuler(euler_2_front'.*degree,cs);
            ori_2_back = orientation.byEuler(euler_2_back'.*degree,cs);
            ori_1 = mean([ori_1_front;ori_1_back],cs);         
            ori_2 = mean([ori_2_front;ori_2_back],cs);
            GBC{idx,3} = ori_1;
            GBC{idx,4} = ori_2;
            GBC{idx,5} = angle(ori_1,ori_2,cs)./degree;
            % reconstruct orientation euler angle map
            ori_1_front_mean = mean(ori_1_front,cs);
            ori_1_back_mean = mean(ori_1_back,cs);
            ori_2_front_mean = mean(ori_2_front,cs);
            ori_2_back_mean = mean(ori_2_back,cs);
            EUmap_front(:,grain_1_front) = repmat([ori_1_front_mean.phi1/degree;...
                ori_1_front_mean.Phi/degree; ori_1_front_mean.phi2/degree],[1,sum(grain_1_front,'all')]);
            EUmap_front(:,grain_2_front) = repmat([ori_2_front_mean.phi1/degree;...
                ori_2_front_mean.Phi/degree; ori_2_front_mean.phi2/degree],[1,sum(grain_2_front,'all')]);
            EUmap_back(:,grain_1_back) = repmat([ori_1_back_mean.phi1/degree;...
                ori_1_back_mean.Phi/degree; ori_1_back_mean.phi2/degree],[1,sum(grain_1_back,'all')]);
            EUmap_back(:,grain_2_back) = repmat([ori_2_back_mean.phi1/degree;...
                ori_2_back_mean.Phi/degree; ori_2_back_mean.phi2/degree],[1,sum(grain_2_back,'all')]);
        end
        workbar(idx/num_gb,sprintf('processing %d / %d boundaries...',[idx num_gb]));
    end
    % EUmap_front = permute(EUmap_front,[2,3,1]);
    % EUmap_back = permute(EUmap_back,[2,3,1]);
end

function [xx, yy] = interp_GB_trace(x, y, num, options)
    % interpolate GB traces from position (x,y) to given number of
    % datapoints num.
    arguments
        x (:,1) double
        y (:,1) double
        num (1,1) double
        options.smooth (1,1) logical = false
    end
    if max(max(abs(diff(x))),max(abs(diff(y)))) > 10
        try
            [x, y, ~] = sort_trace(x, y);
        catch
            warning('GB trace is not ideal')
        end
    end
    num_prev = length(x);
    if num_prev >= num
        xx = x;
        yy = y;
    else
%         idx = linspace(1,num_prev,num);
%         xx = (ceil(idx)-idx)'.*x(floor(idx)) + (idx-floor(idx))'.*x(ceil(idx)); % interpolate every idx
%         xx(idx==ceil(idx)) = x(idx(idx==ceil(idx))); % modify errorous ones with integer idx
%         yy = (ceil(idx)-idx)'.*y(floor(idx)) + (idx-floor(idx))'.*y(ceil(idx));
%         yy(idx==ceil(idx)) = y(idx(idx==ceil(idx)));
        if ~options.smooth
            idx = linspace(1,num_prev,num);
            flag = (idx-floor(idx) <= ceil(idx)-idx)';
            xx = flag.*x(floor(idx)) + (~flag).*x(ceil(idx));
            yy = flag.*y(floor(idx)) + (~flag).*y(ceil(idx));
        else
            idx = linspace(1,num_prev,num);
            frac = idx-floor(idx);
            xx = (1-frac)'.*x(floor(idx)) + frac'.*x(min(num_prev,floor(idx)+1));
            yy = (1-frac)'.*y(floor(idx)) + frac'.*y(min(num_prev,floor(idx)+1));
        end

        
    end
end

function [xx, yy, sorted_list] = sort_trace(x, y)
    % find adjacent relationship
    n = length(x);
    adjPair = {};
    idx = 1;
    for ii = 1:n
        for jj = 1:n
            if ii == jj
                continue
            else
                pos1 = [x(ii) y(ii)];
                pos2 = [x(jj) y(jj)];
                distance = (pos1-pos2)*(pos1-pos2)';
                if distance <= 2
                    adjPair{idx} = [ii, jj];
                    idx = idx + 1;
                end
            end
        end
    end
    % create quasi hash map for adjacent pair
    adjMatrix = zeros(n,3);
    for ii = 1:n
        adjMatrix(ii,1) = ii;
    end
    for ii = 1:length(adjPair)
        idx = adjPair{ii}(1);
        if adjMatrix(idx,2) ~= 0
            adjMatrix(idx,3) = adjPair{ii}(2);
        else
            adjMatrix(idx,2) = adjPair{ii}(2);
        end
    end
    idx = find(adjMatrix(:,3) == 0);
    sorted_list = zeros(1,n);
    num_seg = length(idx) / 2;
    if num_seg == 2
        sorted_list = zeros(2,n);
        for ii = 1:num_seg
            sorted_list(ii,1) = idx(ii*2-1);
            idx = 2;
            while adjMatrix(sorted_list(idx),3)
                if find(sorted_list==adjMatrix(sorted_list(idx-1),2))
                    sorted_list(ii,idx) = adjMatrix(sorted_list(idx-1),3);
                else
                    sorted_list(ii,idx) = adjMatrix(sorted_list(ii-1),2);
                end
                idx = idx + 1;
            end
        end
        sorted_list = linspace(1,n,n);
        % need some more work on this
    elseif num_seg > 2
        sorted_list = linspace(1,n,n);
    else
        sorted_list(1) = idx(1);
        for ii = 2:n
            if find(sorted_list==adjMatrix(sorted_list(ii-1),2))
                sorted_list(ii) = adjMatrix(sorted_list(ii-1),3);
            else
                sorted_list(ii) = adjMatrix(sorted_list(ii-1),2);
            end
        end
    end
    xx = x(sorted_list);
    yy = y(sorted_list);
end

function plot_boundary(boundary, idx)
    % this is the function to plot boundary no.idx together with two
    % adjacent grains
    if idx < 0 || idx > boundary.NumOfBoundaries
        error("boundary number out of range")
    end
    grain_idx = boundary.grainAdjacent{idx};
    map = zeros(size(boundary.grainIdxMap));
    map(boundary.grainIdxMap==grain_idx(1)) = 1;
    map(boundary.grainIdxMap==grain_idx(2)) = 2;
    map(boundary.PixelIdxList{idx}) = 3;
    [row, col] = find(map>0);
    xmin = min(row);
    xmax = max(row);
    ymin = min(col);
    ymax = max(col);
    figure, imshow(map,[0 3],'border','tight')
    figure, imshow(imcrop(map,[ymin xmin ymax-ymin xmax-xmin]),[0 3], 'border','tight')
    colormap(jet)
end

function [xticks, yticks, xticklabels, yticklabels] = generate_tick(x_len, map_size, options)
    arguments
        x_len (1,1) double
        map_size (1,2) double
        options.interval (1,1) double = 10 % mm
    end
    x_pxl = map_size(1);
    y_pxl = map_size(2);
    % length for each pixel
    pxl_len = x_len / x_pxl; % in unit of mm / pxl
    y_len = x_len * y_pxl / x_pxl;
    x_interval = fix(x_len / options.interval);
    y_interval = fix(y_len / options.interval);
    xticks = linspace(0,options.interval*x_interval/pxl_len,x_interval+1);
    xticklabels = num2str(linspace(0,options.interval*x_interval,x_interval+1));
    yticks = linspace(0,options.interval*y_interval/pxl_len,y_interval+1);
    yticklabels = num2str(linspace(0,options.interval*y_interval,y_interval+1));
end

% NOT FINISHED
function [plane] = trace2plane(trace_top, trace_bot, thickness)
% the function is to covert two-layer traces into three-dimensional plane
% info with provided thickness
    arguments
        trace_top (:,2) double
        trace_bot (:,2) double
        thickness (1,1) double = 2
    end
    % preprocess input trace data to be the same length
    len_top = length(trace_top);
    len_bot = length(trace_bot);
    if len_top < len_bot
        trace_top = interp_GB_trace(trace_top(:,1),trace_top(:,2),len_bot);
    elseif len_top > len_bot
        trace_bot = interp_GB_trace(trace_bot(:,1),trace_bot(:,2),len_top);
    end
    % need a step to identify the start-finish order for a trace
    plane = zeros(length(trace_top), thickness, 3);
    for ii = 1:thickness
        for jj = 1:length(trace_top)
            x = trace_bot(jj,1);
            y = trace_bot(jj,2);
%             plane(jj,ii,:) = [trace_bot(, ,ii]
        end
    end
end