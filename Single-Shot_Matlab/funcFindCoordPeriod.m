function [coord00 coord01 coord10 periodHorizontal periodVertical]=funcFindCoordPeriod(image,wildGuess)

global rows columns 


imageFFT=ifft2(image);%FFT(fast fourier transform) in matlab


imageFFT=circshift(imageFFT,[0 round(columns/2)]);
imageFFT=transpose(circshift(transpose(imageFFT),[0 round(rows/2)]));%imageFFT=fftshift(imageFFT);
rowMin=round(rows/2)-10;rowMax=round(rows/2)+10;
columnMin=round(columns/2)-10;columnMax=round(columns/2)+10;
sub00=imageFFT(rowMin:rowMax,columnMin:columnMax);
abssub00=abs(sub00);
peakAmplitude00=max(abssub00(:));
[i,j]=find(abs(sub00)==peakAmplitude00);
coord00=[i,j]+[rowMin,columnMin]-[1,1];
columnMin=round(columns/2)+wildGuess;columnMax=columns;
sub10=imageFFT(coord00(1,1),columnMin:columnMax);
abssub10=abs(sub10);
peakAmplitude10=max(abssub10(:));
coord10=[coord00(1,1) 0];
[m,n]=find(abs(sub10)==peakAmplitude10);
coord10(1,2)=n+columnMin-1;
rowMin=1;rowMax=round(rows/2)-wildGuess;
sub01=imageFFT(rowMin:rowMax,coord00(1,2));
abssub01=abs(sub01);
peakAmplitude01=max(abssub01(:));
coord01=[0 coord00(1,2)];
[p,q]=find(abs(sub01)==peakAmplitude01);
coord01(1,1)=p+rowMin-1;
periodHorizontal=norm(coord00-coord10);
periodVertical=norm(coord00-coord01);



