function [g1 g2]=funcPlotImageLineProbe(image,dimension)
subplot(1,2,1)
g1=imresize(imadjust(imagesc(image)),256);
colormap(gray)

subplot(1,2,2)
g2=plot(image(round(dimension/2),:));

