function sinogramout=funcStreakRemovalFFT(sinogram,bandWidth)
% [numberRays,numberAngles,streaks,streaksHighPassFilter,output]

global NZ

[numberRays,numberAngles]=size(sinogram);
streaks=sum(sinogram,2)/numberAngles;
streaksHighPassFilter=imfilter(streaks,bandWidth);%use image filter command to show 'high pass filter' in Mathematica; conv2?

sinogramout=zeros(size(sinogram));
for i=1:NZ
sinogramout(:,i)=sinogram(:,i)-streaksHighPassFilter;
end