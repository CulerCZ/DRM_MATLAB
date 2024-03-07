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
%% background subtratction applied
% convert into drps' stack
% drp_original is in size of [n1 x n2 x th_num x ph_num]
drp_original_top = igray2drp(igray_sample_top,phitheta,exp_para);
drp_original_bot = igray2drp(igray_sample_bot,phitheta,exp_para);
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
exportgraphics(gcf,"D:\DRM\rawData\groupB_Sample04_top\ipf_top.tif",Resolution=300)
close all
figure, imshow(plot_ipf_map(index_result_bot.EUmap),'Border','tight')
exportgraphics(gcf,"D:\DRM\rawData\groupB_Sample04_top\ipf_bot.tif",Resolution=300)
close all

%% filter non-indexed pixels
% background pixels
[n1,n2] = size(drp_original);
% non_index = ones(n1,n2,'logical'); % 1 for indexed; 0 for non-indexed
index_sum = zeros(n1,n2);
for ii = 1:n1
    for jj = 1:n2
        index_num(ii,jj) = sum(drp_original{ii,jj},'all');
    end
end
non_index_bg = index_num > 3e4;
figure, imshow(non_index_bg)
% poorly-indexed pixels
non_index_poor = index_result.quality < prctile(index_result.quality,95);

non_index = non_index_poor & non_index_bg; 
figure, imshow(non_index)
%% plot results
% [ebsd_DRM] = plot_ipf_DRM(index_result.EUmap, 'Ni', nonindex=ones(size(index_result.EUmap)));
figure, imshow(plot_ipf_map(index_result_bot.EUmap))
%% check indexing results
[drp_measurement, drp_predicted, x, y] = check_indexing_result(index_result.EUmap,drp_original,exp_para);
%%
drp_out = drp_encode(AE_DRM,drp_measurement,exp_para);
drp_tmp = drp_encode(AE_DRM,drp_predicted,exp_para);

%% direct indexing without autoencoder
bandIntensity = [0.1 0.3 0.5 0.7];
for ii = 1:4
    exp_para.fitting_para(2) = bandIntensity(ii);
    drpLib = DRPLibGenerator(5*degree, exp_para);
    
    indexResult = DirectDIEngine(drp_original, drpLib);
    figure, imshow(plot_ipf_map(indexResult.euMap));
    title(sprintf("band peak ratio %.1f",bandIntensity(ii)),'FontSize',14,'FontWeight','bold')
end
%% generate DRP dictionary
drpLib = DRPLibGenerator(3*degree, exp_para);
%% dictionary indexing
% indexResult_top = DirectDIEngine(drp_original_top, drpLib);
indexResult_bot = DirectDIEngine(drp_original_bot,drpLib);
%% check indexing outcome
check_indexing_result(indexResult_top.euMap,drp_original_top,exp_para);


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