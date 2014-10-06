function [g1 g2]=funcPlotArrayLineProbe(image,dimension,zMin,zMax)
subplot(1,2,1)
%g1=imagesc(zMin,zMax,image);
clims=[zMin zMax];
g1=imresize(imadjust(imagesc(image,clims)),256);
colormap(gray)

subplot(1,2,2)
g2=plot(image(round(dimension/2),:));