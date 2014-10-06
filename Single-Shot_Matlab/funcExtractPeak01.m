function [peak01,rows01,columns01]=funcExtractPeak01(image)

global coord01 periodHorizontal periodVertical
%coord01=[513 886];
%periodHorizontal=245;
%periodVertical=196;

rowMin=coord01(1,1)-floor(periodHorizontal/2);
rowMax=coord01(1,1)+floor(periodHorizontal/2);
columnMin=coord01(1,2)-floor(periodVertical/2);
columnMax=coord01(1,2)+floor(periodVertical/2);
peak01=image(rowMin:rowMax,columnMin:columnMax);
[rows01,columns01]=size(peak01);
