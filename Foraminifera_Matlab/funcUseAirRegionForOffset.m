function [sinogramout]=funcUseAirRegionForOffset(sinogram,airRegionWidth)
%[numberRays,numberAngles,airTop,airBottom,offset]
[numberRays,numberAngles]=size(sinogram);
airTop=sinogram(1:airRegionWidth,:);
airTop=sum(airTop,1)/airRegionWidth;
airTop=transpose(airTop);
airBottom=sinogram(numberRays-airRegionWidth+1:numberRays,:);
airBottom=sum(airBottom,1)/airRegionWidth;
airBottom=transpose(airBottom);
offset=(airTop+airBottom)/2;
sinogramtrans=transpose(sinogram);
sinogramout=zeros(size(sinogramtrans));
for i=1:numberRays
sinogramout(:,i)=transpose(sinogramtrans(:,i)-offset);
end
sinogramout=transpose(sinogramout);