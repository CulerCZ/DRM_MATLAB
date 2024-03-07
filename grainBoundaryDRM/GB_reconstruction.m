% figure, hold on
thickpxlratio = ceil(size(drp_original_top,1)/15);
GB_voxel = cell(length(GBC),3);
for grainNum = 1:length(GBC)
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

    [gbx, gby, gbz] = ind2sub(size(voxels),find(voxels<=1));
    gbVoxNum = length(gbx);
    start_voxel = [xmin, ymin, 0];
    % for ii = 1:gbVoxNum
    %     % plotCube(start_voxel+[gbx(ii),gby(ii),gbz(ii)],1,color=[173, 216, 230]/256,facealpha=0.717)
    %     % scatter3(start_voxel(1)+gbx(ii),)
    % end
    % scatter3(start_voxel(1)+gbx, start_voxel(2)+gby, start_voxel(3)+gbz, 2, "blue")
    GB_voxel{grainNum,1} = start_voxel(1)+gbx;
    GB_voxel{grainNum,2} = start_voxel(2)+gby;
    GB_voxel{grainNum,3} = start_voxel(3)+gbz;
    fprintf("grain boundary %d finished.\n",grainNum)
end

%% plot GB network in 3d
cmap = colormap(jet);

figure, hold on
for ii = 1:length(GBC)
    color_rgb = cmap(fix(GBC{ii,5}/63*256),:);
    scatter3(GB_voxel{ii,1},GB_voxel{ii,2},GB_voxel{ii,3},1,color_rgb)

    % blue GB traces
    % scatter3(GB_voxel{ii,1},GB_voxel{ii,2},GB_voxel{ii,3},1,"blue")

end

axis equal
xlim([-1 size(drp_original_top,1)+1]) 
zlim([-1 thickpxlratio+1])
ylim([-1 size(drp_original_top,2)+1])
set(gca,'visible','off')


img_front = plot_ipf_map(index_result_top.EUmap);
img_back = plot_ipf_map(index_result_bot.EUmap);
[width, height, ~] = size(img_front);

zPlane = 0;
x = [0 width;0 width];
y = [0 0; height height];
z = [zPlane,zPlane;zPlane,zPlane];
img_back = cat(3,img_back(:,:,1).',img_back(:,:,2).',img_back(:,:,3).');
surf(x,y,z,img_back,'FaceColor','texturemap',FaceAlpha=0.3,EdgeColor='none')

zPlane = thickpxlratio;
x = [0 width;0 width];
y = [0 0; height height];
z = [zPlane,zPlane;zPlane,zPlane];
img_front = cat(3,img_front(:,:,1).',img_front(:,:,2).',img_front(:,:,3).');

surf(x,y,z,img_front,'FaceColor','texturemap',FaceAlpha=0.3,EdgeColor='none')

%%
function plotCube(center, l, options)
    arguments
        center (1,3) double
        l (1,1) double
        options.color (1,3) double = [77 77 77]/256
        options.facealpha (1,1) double = 0.6
    end
    V1=[-1;  1; 1; -1; -1;  1; 1; -1;];
    V2=[-1; -1; 1;  1; -1; -1; 1;  1;];
    V3=[-1; -1;-1; -1;  1;  1; 1;  1;];
    F=[1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8;];
    Vertices=[V1 V2 V3]*l/2 + center;
    patch(gca,'Faces',F,'Vertices',Vertices,'FaceColor',options.color, ...
        'FaceAlpha',options.facealpha,'EdgeColor','k','LineWidth',1); 

end