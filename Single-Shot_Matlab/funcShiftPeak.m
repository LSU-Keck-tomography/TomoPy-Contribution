function [shiftedPeak]=funcShiftPeak(peak,rows,columns)
shiftedPeak=zeros(rows,columns);
[rowsPeak,columnsPeak]=size(peak);
rShift=round((rows-rowsPeak)/2);
cShift=round((columns-columnsPeak)/2);
for r=1:rowsPeak
    for c=1:columnsPeak
shiftedPeak(r+rShift,c+cShift)=peak(r,c);%Table form
    end
end