cs = crystalSymmetry('cubic');
ori_1deg = equispacedSO3Grid(cs,'resolution',1*degree);

% in a grain with orientation, the direction in crystal reference frame has
% the corresponding direction in specimen reference system
r_001 = ori_1deg * Miller({0,0,1},cs);
r_100 = ori_1deg * Miller({1,0,0},cs);
r_010 = ori_1deg * Miller({0,1,0},cs);

figure,
scatter3(r_001.x, r_001.y, r_001.z, 5,'red','filled')
hold on
% scatter3(r_010.x, r_010.y, r_010.z, 5,'blue','filled')
% scatter3(r_100.x, r_100.y, r_100.z, 5,'green','filled')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'-','LineWidth',3,'Color','black')
axis equal
xlim([-1.3 1.3])
ylim([-1.3 1.3])
zlim([-1.3 1.3])
set(gca,'Visible','off')
text(1.3,0,0,'X',"FontSize",14)
text(0,1.3,0,'Y',"FontSize",14)
text(0,0,1.3,'Z',"FontSize",14)
view(130,15)

% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/peakFZ_001.tif",'Resolution',600)

%% plot peak positions
th_fz = r_001.theta;
ph_fz = r_001.rho;
th_fz_r = 2*th_fz;
x_fz = sin(th_fz_r).*cos(ph_fz);
y_fz = sin(th_fz_r).*sin(ph_fz);
z_fz = cos(th_fz_r);
pos_xyz = unique([x_fz, y_fz, z_fz],'row');
figure, hold on
scatter3(pos_xyz(:,1),pos_xyz(:,2),pos_xyz(:,3),5,'black','filled')
scatter3(r_001.x*0.5, r_001.y*0.5, r_001.z*0.5, 5,'red','filled')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'-','LineWidth',3,'Color','black')
axis equal
xlim([-1.3 1.3])
ylim([-1.3 1.3])
zlim([-1.3 1.3])
set(gca,'Visible','off')
text(1.3,0,0,'X',"FontSize",14)
text(0,1.3,0,'Y',"FontSize",14)
text(0,0,1.3,'Z',"FontSize",14)
view(130,15)
exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/refPeakFZ_001.tif",'Resolution',600)


%% plot reflectance peak position for all three
th_fz = r_001.theta;
ph_fz = r_001.rho;
th_fz_r = 2*th_fz;
x_fz = sin(th_fz_r).*cos(ph_fz);
y_fz = sin(th_fz_r).*sin(ph_fz);
z_fz = cos(th_fz_r);
pos_xyz = unique([x_fz, y_fz, z_fz],'row');
figure, hold on
scatter3(pos_xyz(:,1),pos_xyz(:,2),pos_xyz(:,3),5,'red','filled')

th_fz = r_100.theta;
ph_fz = r_100.rho;
th_fz_r = 2*th_fz;
valid = th_fz<=pi/2;
x_fz = sin(th_fz_r).*cos(ph_fz);
y_fz = sin(th_fz_r).*sin(ph_fz);
z_fz = cos(th_fz_r);
pos_xyz = unique([x_fz(valid), y_fz(valid), z_fz(valid)],'row');
scatter3(pos_xyz(:,1),pos_xyz(:,2),pos_xyz(:,3),5,'green','filled')

th_fz = r_010.theta;
ph_fz = r_010.rho;
th_fz_r = 2*th_fz;
valid = th_fz<=pi/2;
x_fz = sin(th_fz_r).*cos(ph_fz);
y_fz = sin(th_fz_r).*sin(ph_fz);
z_fz = cos(th_fz_r);
pos_xyz = unique([x_fz(valid), y_fz(valid), z_fz(valid)],'row');
scatter3(pos_xyz(:,1),pos_xyz(:,2),pos_xyz(:,3),5,'blue','filled')

quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'-','LineWidth',3,'Color','black')
axis equal
xlim([-1.3 1.3])
ylim([-1.3 1.3])
zlim([-1.3 1.3])
set(gca,'Visible','off')
text(1.3,0,0,'X',"FontSize",14)
text(0,1.3,0,'Y',"FontSize",14)
text(0,0,1.3,'Z',"FontSize",14)
view(130,15)
% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/refPeakFZ_all.tif",'Resolution',600)


cs = crystalSymmetry('cubic');
ori_1deg = equispacedSO3Grid(cs,'resolution',2*degree);

% in a grain with orientation, the direction in crystal reference frame has
% the corresponding direction in specimen reference system
r_111 = ori_1deg * Miller({1,1,1},cs);
r_100 = ori_1deg * Miller({-1,1,1},cs);
r_010 = ori_1deg * Miller({1,-1,1},cs);
r_000 = ori_1deg * Miller({-1,-1,1},cs);
xyz111 = normr(r_111.xyz);
xyz100 = normr(r_100.xyz);
xyz010 = normr(r_010.xyz);
xyz000 = normr(r_000.xyz);
figure,
scatter3(xyz111(:,1), xyz111(:,2), xyz111(:,3), 5,'red','filled')
hold on
scatter3(xyz100(:,1), xyz100(:,2), xyz100(:,3), 5,'blue','filled')
scatter3(xyz010(:,1), xyz010(:,2), xyz010(:,3), 5,'green','filled')
scatter3(xyz000(:,1), xyz000(:,2), xyz000(:,3), 5,'yellow','filled')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'-','LineWidth',3,'Color','black')
axis equal
xlim([-1.3 1.3])
ylim([-1.3 1.3])
zlim([-1.3 1.3])
set(gca,'Visible','off')
text(1.3,0,0,'X',"FontSize",14)
text(0,1.3,0,'Y',"FontSize",14)
text(0,0,1.3,'Z',"FontSize",14)
view(125,15)

% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/peakFZ_111only.tif",'Resolution',600)

%% reflectance peak distribution
th_fz = r_111.theta;
ph_fz = r_111.rho;
th_fz_r = 2*th_fz;
x_fz = sin(th_fz_r).*cos(ph_fz);
y_fz = sin(th_fz_r).*sin(ph_fz);
z_fz = cos(th_fz_r);
pos_xyz = unique([x_fz, y_fz, z_fz],'row');
figure, hold on
scatter3(pos_xyz(:,1),pos_xyz(:,2),pos_xyz(:,3),5,'black','filled')
scatter3(xyz111(:,1)*0.5, xyz111(:,2)*0.5, xyz111(:,3)*0.5, 5,'red','filled')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'-','LineWidth',3,'Color','black')
axis equal
xlim([-1.3 1.3])
ylim([-1.3 1.3])
zlim([-1.3 1.3])
set(gca,'Visible','off')
text(1.3,0,0,'X',"FontSize",14)
text(0,1.3,0,'Y',"FontSize",14)
text(0,0,1.3,'Z',"FontSize",14)
view(40,15)
exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/refPeakFZ_111.tif",'Resolution',600)
