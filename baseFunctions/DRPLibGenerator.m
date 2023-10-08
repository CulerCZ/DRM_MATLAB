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