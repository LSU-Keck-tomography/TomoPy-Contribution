function [peak10,rows10,columns10]=funcExtractPeak10(image)

global coord10 periodHorizontal periodVertical
%coord10=[513 886];
%periodHorizontal=245;
%periodVertical=196;

rowMin=coord10(1,1)-floor(periodHorizontal/2);
rowMax=coord10(1,1)+floor(periodHorizontal/2);
columnMin=coord10(1,2)-floor(periodVertical/2);
columnMax=coord10(1,2)+floor(periodVertical/2);
peak10=image(rowMin:rowMax,columnMin:columnMax);
[rows10,columns10]=size(peak10);
