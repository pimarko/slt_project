workingDir = '/Users/Nico/Documents/Matlab/SLT/pyhon/gif';
outputVideo = VideoWriter(fullfile(workingDir,'video.avi'));
outputVideo.FrameRate = 2;
open(outputVideo)

for ii = 2:73
   img = imread(fullfile(workingDir,sprintf('plot_final_clusters_%d.png', ii)));
   writeVideo(outputVideo,img)
end

close(outputVideo)