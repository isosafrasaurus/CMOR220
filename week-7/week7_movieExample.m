% Anastasiya Protasov, CMOR220, Fall 2023, movie example
% week7_movieExample.m
% Script to show and example of creating an animation in MATLAB (lorenz attractor)
% list of input: none
% list of output: none
% Last Modified: OCtober 4, 2023

%the following line creates the video file, you will need to modify this 
%in order to choose a path on your own machine on which to save your video 
%myVid = VideoWriter('/Users/anastasiyaprotasov/Desktop/CAAM_CAAM210/Spring2023/MATLAB/Week6.2,''.avi');
% myVid = VideoWriter('Example.avi');
% open(myVid)     % open your video file
% 
% for j = 1:length(t2)
%     %plot the entire solution as well as a single point of the solution 
%     plot3(y2(:,1),y2(:,2),y2(:,3),y2(j,1),y2(j,2),y2(j,3),'ro');  
%     xlabel('x');ylabel('y');zlabel('z');
%     title(['lorenz attractor, time = ', num2str(t2(j))]);
%     currFrame = getframe(gcf);         % save your figure as a frame
%     writeVideo(myVid,currFrame);       % add the frame to your movie
% end
% 
% close(myVid);   % close the video file, This must be done for the video
%                 %file to play properly


% make another video


otherVid = VideoWriter('Example.avi');
open(otherVid);

for j = 1:length(t2)
    plot(y2(:,1), y2(:,3), y2(j,1),y2(j,3),'ro');  %projection on the xz-plane
    xlabel('x'); ylabel('z')
    title(['lorenz attractor on the xz-plane']);
    currFrame = getframe(gcf);
    writeVideo(otherVid,currFrame);
end

close(otherVid);
    
%make another video

% thirdVid = VideoWriter('Example.avi');
% open(thirdVid);
% 
% for j = 1:length(t2)
%     plot(y2(1:j,1), y2(1:j,3), y2(j,1),y2(j,3),'ro');  % draw as you go
%     axis([-20,20,0,50])
%     xlabel('x'); ylabel('z');
%     title(['lorenz attractor on the xz-plane']);
%     currFrame = getframe(gcf);
%     writeVideo(thirdVid,currFrame);
% end
% 
% close(thirdVid);
