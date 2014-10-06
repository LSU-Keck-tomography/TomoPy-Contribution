function [peak00,rows00,columns00]=funcExtractPeak00(image)

global coord00 periodHorizontal periodVertical
%coord00=[513 886];
%periodHorizontal=245;
%periodVertical=196;

rowMin=coord00(1,1)-floor(periodHorizontal/2);
rowMax=coord00(1,1)+floor(periodHorizontal/2);
columnMin=coord00(1,2)-floor(periodVertical/2);
columnMax=coord00(1,2)+floor(periodVertical/2);
peak00=image(rowMin:rowMax,columnMin:columnMax);
[rows00,columns00]=size(peak00);
