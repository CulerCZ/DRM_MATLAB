%% 3d animation 
filename = 'GB_network_animation.gif';
nFrames = 36;
angleStep = 360/nFrames;

for ii = 1:nFrames
    view(angleStep*ii, 30);
    drawnow

    frame = getframe(gcf);
    img = frame2im(frame);
    [A,map] = rgb2ind(img,256);

    if ii == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1)
    else
        imwrite(A,map,filename,'WriteMode','append','DelayTime',0.1)
    end
end
