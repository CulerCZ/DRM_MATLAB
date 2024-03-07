% Parameters
radius = 0.05; % Radius of the cylinder
length = 0.8; % "Height" of the cylinder, but it's along the x-axis in this orientation

% Angles for the semicircle (half cylinder)
theta = linspace(0, pi, 50);

% Length along the cylinder's axis (x-axis in this case)
X = linspace(-length/2, length/2, 50);
 
% Generate mesh for cylinder
[Theta, X] = meshgrid(theta, X); % Create a grid of theta and x values
Y = radius * cos(Theta); % Y coordinates
Z = radius * sin(Theta); % Z coordinates, upwards
C = Z./Y;
% Plot the half cylinder
figure;
surf(X, Y, Z, 'LineStyle', 'none'); % Plotting the surface
colormap("gray")
hold on; % Keep the plot for adding more graphics
axis equal; % Keep the aspect ratio equal
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Half Cylinder on xOy Plane with Height Along x-axis');

% Define and plot the rectangle (base)
rectX = [-length/2, length/2, length/2, -length/2];
rectY = [-radius, -radius, radius, radius];
rectZ = [0, 0, 0, 0]; % Since it's on the xOy plane, Z is 0 for all points
fill3(rectX, rectY, rectZ, 'k', 'FaceAlpha', 0.5); % Adding the rectangle as the base

% Fill the half-circles at each end
theta_fill = linspace(0, pi, 100); % More points for a smoother half-circle
x0 = -length/2; % Starting x position for the first half-circle
x1 = length/2; % Ending x position for the second half-circle

% Coordinates for the first half-circle (start)
y0 = radius * cos(theta_fill);
z0 = radius * sin(theta_fill);
fill3(x0*ones(size(theta_fill)), y0, z0, 'k', 'FaceAlpha', 0.5);

% Coordinates for the second half-circle (end)
y1 = radius * cos(theta_fill);
z1 = radius * sin(theta_fill);
fill3(x1*ones(size(theta_fill)), y1, z1, 'k', 'FaceAlpha', 0.5);

% Parameters for the half circle
radius = 0.8;
theta = linspace(0, pi, 100); % Half circle from 0 to pi

% Coordinates for the half circle on yOz plane
Y = radius * cos(theta); % Y-axis
Z = radius * sin(theta); % Z-axis
y1 = radius * cos(theta_fill);
z1 = radius * sin(theta_fill);
fill3(zeros(size(theta_fill)), y1, z1, 'r', 'FaceAlpha', 0.5);
% Plotting
plot3(zeros(size(Y)), Y, Z, 'r', 'LineWidth',2); % Plotting in blue

% Plot 3d dRP
exp_para_plot.th_max = 65;
exp_para_plot.th_min = 8;
exp_para_plot.th_num = 58;
exp_para_plot.ph_min = 0;
exp_para_plot.ph_max = 359;
exp_para_plot.ph_num = 360;
exp_para_plot.faceting = [1 0 0];
exp_para_plot.fitting_para = [1, 0.6, 18, 6, 0.8, 8];
drpsim_tmp = DRPsim(-5,45,0,exp_para_plot); % the matrix in form of th_num * ph_num
theta = repmat(linspace(exp_para_plot.th_min, exp_para_plot.th_max, exp_para_plot.th_num)',1,exp_para_plot.ph_num+1);
phi = repmat(0:360/exp_para_plot.ph_num:360,exp_para_plot.th_num,1);
x = cosd(theta).*cosd(phi);
y = cosd(theta).*sind(phi);
z = sind(theta);
CC = [drpsim_tmp, drpsim_tmp(:,1)];
surf(x,y,z,CC,'EdgeColor','none','FaceAlpha',0.4);
axis equal
colormap('jet')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'k-','LineWidth',2,'MaxHeadSize',0.6)
% zlim([zsurf 1.3])
% xlim([-1 1.3])
% ylim([-1 1.3])
view(130,20)

hold off; % Release the plot
set(gca,"Visible",'off')
exportgraphics(gcf,"/Users/chenyangzhu/Desktop/curveDRMfigures/edge3Dschematic.tif",Resolution=300)