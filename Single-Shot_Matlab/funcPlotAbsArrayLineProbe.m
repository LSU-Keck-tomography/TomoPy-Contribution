function [g1 g2]=funcPlotAbsArrayLineProbe(image,dimension)
absImage=abs(image);
subplot(1,2,1)
g1=imresize(imadjust(imagesc(real(absImage))),256);%Note: image is a complex matrix
colormap(gray)

subplot(1,2,2)
g2=plot(absImage(round(dimension/2),:));
