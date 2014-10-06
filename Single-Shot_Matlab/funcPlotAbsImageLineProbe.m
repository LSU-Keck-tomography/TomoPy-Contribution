function [g1 g2]=funcPlotAbsImageLineProbe(image,dimension)
absImage=abs(image);
subplot(1,2,1)
g1=imresize(imadjust(imagesc(real(image))),256);%Note: image is a complex matrix
colormap(gray)

subplot(1,2,2)
g2=plot(absImage(round(dimension/2),:));
colormap(gray)

