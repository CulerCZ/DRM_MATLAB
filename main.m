%% DRM measurement indexing engine
exp_para.th_max = 65;
exp_para.th_min = 8;
exp_para.th_num = 20;
exp_para.ph_num = 72;
exp_para.ph_min = 0;
exp_para.ph_max = 355;
exp_para.faceting = [1 0 0];
exp_para.fitting_para = [1, 0.6, 15, 8, 0.8, 8];

%% load sample and background dataset
[igray_sample, phitheta, pos, img_sample] = drp_loader(exp_para,zeros(1,4),format='jpg',scale=0.2);
[igray_back, ~, ~] = drp_loader(exp_para,pos,scale=0.2);

%% background subtratction applied
igray_norm = bg_subtraction(igray_sample,igray_back,1.5);
[n1,n2,n3] = size(igray_norm);

% convert into drps' stack
% drp_original is in size of [n1 x n2 x th_num x ph_num]
drp_original_02 = igray2drp(igray_norm,phitheta,exp_para);
% uisave('igray_norm','igray_norm')
%% clear igray_norm igray_sample igray_back
drp_sample = igray2drp(igray_sample,phitheta,exp_para);

%% show DRP
% x = 550;
% y = 1500;
% figure, DRPdisp(drp_original{600,600},exp_para);
% clear x y
drp_measurement = check_measurement(img_sample,drp_original_02,exp_para); 
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
tic
% drp_rescale = cellfun(@(x) uint8(rescale(x)*255), drp_sample, 'UniformOutput', false);
index_result = IndexingEngine(drp_original,AE_DRM,exp_para,drpDic,euDic,rotDic);
toc

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
figure, imshow(plot_ipf_map(index_result.EUmap))
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
%%
% drpLib_sparce = DRPLibGenerator(5*degree, exp_para);

indexResult_layer01 = DirectDIEngine(drp_original_02, drpLib);


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