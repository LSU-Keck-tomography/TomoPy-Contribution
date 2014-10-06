function sinogram=funcShiftSinogram(sinogram,shift,numberAngles,numberColumns)
%[rInitial, cInitial, rFinal, cFinal, sinoShifted,error]
[rInitial,cInitial]=size(sinogram);
transposeSinogram=transpose(sinogram);
if shift>0
    sinoShifted=transposeSinogram(:,(shift:numberColumns));
    [r,c]=size(sinoShifted);
    sinoShifted=padarray(sinoShifted,[numberAngles-r numberColumns-c],'post');
    sinoShifted=transpose(sinoShifted);
elseif shift<0
    sinoShifted=transposeSinogram(:,(1:numberColumns+shift));
    [r,c]=size(sinoShifted);
    sinoShifted=padarray(sinoShifted,[numberAngles-r numberColumns-c],'post');
    sinoShifted=transpose(sinoShifted);
else
    sinoShifted=sinogram;
end

[rFinal,cFinal]=size(sinoShifted);
if rInitial==rFinal && cInitial==cFinal
    error=0;
else error=1;
end
sinogram=sinoShifted;