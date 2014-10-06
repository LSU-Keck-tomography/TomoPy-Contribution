%Exp17_LatexSpheres_2DGr_ver1.m
%Single shot phase contrast with a single 2D grating. The grating is a checkerboard with Pi-phase shift at 18 keV.
%This program makes a lot of figures and this makes the filesize large and increases the time to save the file.  
%To reduce filesize (at the cost of erasing graphics), on the Mathematica menu bar, select Cell/Delete All Output

%% Initialization
clear ; close all; clc
%set path directory
%cd('/Users/jumaoyuan/Desktop/phase_contrast/temp_large_files/Exp17_225um_PS_sphere_2D_18kev_4.8um')

%%
%------------------------------------------------------------------------------------------------
%How many rows,columns in the images?  Filenames for ref, dark, sample files
%------------------------------------------------------------------------------------------------

filenamesReference=dir('*_ref_*.tif');%read reference files 
filenamesDark=dir('*_dark_*.tif');%read dark-field files
filenamesSample=dir('*_sample_*.tif');%read sample files
fprintf('Number of sample files.\n');
length(filenamesSample)%number of sample files
filenamesSample(1).name;%first directory of sample files
fprintf('Program paused. Press enter to continue.\n');
pause;
%filenamesReference.name
%filenamesDark.name
%filenamesSample.name


%%
%-------------------------------------------------------------------------------------------------
%Calculate average image:   dataDark, dataRef  Find coordinates of (1,0), (0,1), and (0,0) peaks.
%-------------------------------------------------------------------------------------------------
dataDark=funcCalculateAverageImageTIFF(filenamesDark);%call Function:funcCalculateAverageImageTIFF
dataRef=funcCalculateAverageImageTIFF(filenamesReference);
[rows,columns]=size(dataRef);

global coord00 coord01 coord10 periodHorizontal periodVertical
[coord00 coord01 coord10 periodHorizontal periodVertical]=funcFindCoordPeriod(dataRef-dataDark,150);
fprintf('values of coord00, coord01, coord10, periodHorizontal, periodVertical ... \n');
PrintCoord=[coord00 coord01 coord10 periodHorizontal periodVertical]
%coord00(1,1)
%coord00(1)
fprintf('Program paused. Press enter to continue.\n');
pause;


image=dataRef-dataDark;
imageFFT=ifft2(image);
imageFFT=circshift(imageFFT,[0 round(columns/2)]);
imageFFT=transpose(circshift(transpose(imageFFT), [0 round(rows/2)]));%imageFFT=fftshift(imageFFT); 


absImageFFT=abs(imageFFT);
figure(1);
subplot(2,2,1)
plot(absImageFFT(coord00(1),:))%coord00(1,1)=513;
subplot(2,2,2)
plot(absImageFFT(coord00(1),:),'*')
subplot(2,2,3)
plot(absImageFFT(:,coord00(2)))%coord00(1,2)=641;
subplot(2,2,4)
plot(absImageFFT(:,coord00(2)),'*')
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%-------------------------------------------------------------------------------------------------
%Process Reference (1,0), (0,1), and (0,0) peaks in Fourier space 
%-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Calculate Reference harmonic image;  filtering is optional.
%--------------------------------------------------------------------------
reference=dataRef-dataDark;
figure(2);
[g1 g2]=funcPlotImageLineProbe(reference,rows);
referenceFFT=ifft2(reference);
referenceFFT=circshift(referenceFFT, [0 round(columns/2)]);
referenceFFT=transpose(circshift(transpose(referenceFFT),[0 round(rows/2)]));%referenceFFT=fftshift(referenceFFT)
%image=referenceFFT is complex 
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%2. Extract (1,0) peak
%--------------------------------------------------------------------------
[peak10,rows10,columns10]=funcExtractPeak10(referenceFFT);
figure(3);
funcPlotAbsImageLineProbe(peak10,rows10);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%3. Sine Bell filter (1,0) peak
%--------------------------------------------------------------------------
figure(4);
subplot(1,2,1)
gBeforeSineBell=plot(abs(peak10(1+round(rows10/2),:)),'r');
hold on
plot(real(peak10(1+round(rows10/2),:)),'b');
hold on
plot(imag(peak10(1+round(rows10/2),:)),'g');
peak10=funcSineBellFilter(peak10);
subplot(1,2,2)
gAfterSineBell=plot(abs(peak10(1+round(rows10/2),:)),'r');
hold on
plot(real(peak10(1+round(rows10/2),:)),'b');
hold on
plot(imag(peak10(1+round(rows10/2),:)),'g');
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%4. Shift (1,0) peak to the center of the image
%--------------------------------------------------------------------------
[shiftedPeak10]=funcShiftPeak(peak10,rows,columns);
figure(5);
[g1 g2]=funcPlotAbsImageLineProbe(shiftedPeak10,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%5. Inverse FFT of (1,0) peak
%--------------------------------------------------------------------------
result10Ref=fft2(shiftedPeak10); %fft2=InverseFourier in Mathematica 
figure(6);
[g1 g2]=funcPlotAbsArrayLineProbe(result10Ref,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%6. Extract (0,0) peak
%--------------------------------------------------------------------------
[peak00,rows00,columns00]=funcExtractPeak00(referenceFFT);
figure(7);
funcPlotAbsImageLineProbe(peak00,rows00);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%7. Sine Bell filter (0,0) peak
%--------------------------------------------------------------------------
figure(8);
subplot(1,2,1)
gBeforeSineBell=plot(abs(peak00(1+round(rows00/2),:)),'r');
hold on
plot(real(peak00(1+round(rows00/2),:)),'b');
hold on
plot(imag(peak00(1+round(rows00/2),:)),'g');
peak00=funcSineBellFilter(peak00);
subplot(1,2,2)
gAfterSineBell=plot(abs(peak00(1+round(rows00/2),:)),'r');
hold on
plot(real(peak00(1+round(rows00/2),:)),'b');
hold on
plot(imag(peak00(1+round(rows00/2),:)),'g');
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%8. Place (0,0) peak in center of image
%--------------------------------------------------------------------------
figure(9);
[shiftedPeak00]=funcShiftPeak(peak00,rows,columns);
funcPlotAbsImageLineProbe(shiftedPeak00,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%9. Inverse FFT of (0,0) peak
%--------------------------------------------------------------------------
result00Ref=fft2(shiftedPeak00); 
figure(10);
[g1 g2]=funcPlotAbsArrayLineProbe(result00Ref,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%10. Extract (0,1) peak
%--------------------------------------------------------------------------
[peak01,rows01,columns01]=funcExtractPeak01(referenceFFT);
figure(11);
funcPlotAbsImageLineProbe(peak01,rows01);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%11. Sine Bell filter (0,1) peak
%--------------------------------------------------------------------------
figure(12);
subplot(1,2,1)
gBeforeSineBell=plot(abs(peak01(1+round(rows01/2),:)),'r');
hold on
plot(real(peak01(1+round(rows01/2),:)),'b');
hold on
plot(imag(peak01(1+round(rows01/2),:)),'g');
peak01=funcSineBellFilter(peak01);
subplot(1,2,2)
gAfterSineBell=plot(abs(peak01(1+round(rows01/2),:)),'r');
hold on
plot(real(peak01(1+round(rows01/2),:)),'b');
hold on
plot(imag(peak01(1+round(rows01/2),:)),'g');
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%12. Place (0,1) peak in center of image
%--------------------------------------------------------------------------
[shiftedPeak01]=funcShiftPeak(peak01,rows,columns);
figure(13);
funcPlotAbsImageLineProbe(shiftedPeak01,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%13. Inverse FFT of (0,1) peak
%--------------------------------------------------------------------------
result01Ref=fft2(shiftedPeak01); %use default parameters
figure(14);
[g1 g2]=funcPlotAbsArrayLineProbe(result01Ref,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;





%%
%Process i-th Sample file
%-------------------------------------------------------------------------------------------------
%Process i-th Sample (1,0), (0,1), and (0,0) peaks in Fourier space   
%-------------------------------------------------------------------------------------------------
index=1;
dataSample=imread(filenamesSample(index).name);
figure(15);
imresize(imadjust(imagesc(real(dataSample))),256);
colormap(gray)
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%1. Calculate Reference harmonic image;  filtering is optional.
%--------------------------------------------------------------------------
sample=double(dataSample)-double(dataDark);
%trans=GaussianFilter[trans,4]; 
figure(16);
[g1 g2]=funcPlotImageLineProbe(sample,rows);
fprintf('Program paused. Press enter to continue.\n');
pause;

sampleFFT=ifft2(sample);
sampleFFT=circshift(sampleFFT,[1 round(columns/2)]);%sampleFFT=fftshift(sampleFFT);
sampleFFT=transpose(circshift(transpose(sampleFFT),[1 round(rows/2)]));
figure(17);
[g1 g2]=funcPlotAbsImageLineProbe(sampleFFT,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%2. Extract (1,0) peak
%--------------------------------------------------------------------------
[peak10,rows10,columns10]=funcExtractPeak10(sampleFFT);
[rows10 columns10];
figure(18);
funcPlotAbsImageLineProbe(peak10,rows10);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%3. Sine Bell filter (1,0) peak
%--------------------------------------------------------------------------
figure(19)
subplot(1,2,1)
gBeforeSineBell=plot(abs(peak10(1+round(rows10/2),:)),'r');
hold on
plot(real(peak10(1+round(rows10/2),:)),'b');
hold on
plot(imag(peak10(1+round(rows10/2),:)),'g');
peak10=funcSineBellFilter(peak10);
subplot(1,2,2)
gAfterSineBell=plot(abs(peak10(1+round(rows10/2),:)),'r');
hold on
plot(real(peak10(1+round(rows10/2),:)),'b');
hold on
plot(imag(peak10(1+round(rows10/2),:)),'g');
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%4. Shift (1,0) peak to the center of the image
figure(20);
[shiftedPeak10]=funcShiftPeak(peak10,rows,columns);
funcPlotAbsImageLineProbe(shiftedPeak10,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%5. Inverse FFT of (1,0) peak
%--------------------------------------------------------------------------
result10Sample=fft2(shiftedPeak10); 
figure(21);
funcPlotAbsArrayLineProbe(result10Sample,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%6. Extract (0,0) peak
%--------------------------------------------------------------------------
[peak00,rows00,columns00]=funcExtractPeak00(sampleFFT);
figure(22);
[g1 g2]=funcPlotAbsImageLineProbe(peak00,rows00);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%7. Sine Bell filter (0,0) peak
%--------------------------------------------------------------------------
figure(23)
subplot(1,2,1)
gBeforeSineBell=plot(abs(peak00(1+round(rows00/2),:)),'r');
hold on
plot(real(peak00(1+round(rows00/2),:)),'b');
hold on
plot(imag(peak00(1+round(rows00/2),:)),'g');
peak00=funcSineBellFilter(peak00);
subplot(1,2,2)
gAfterSineBell=plot(abs(peak00(1+round(rows00/2),:)),'r');
hold on
plot(real(peak00(1+round(rows00/2),:)),'b');
hold on
plot(imag(peak00(1+round(rows00/2),:)),'g');
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%8. Place (0,0) peak in center of image
%--------------------------------------------------------------------------
[shiftedPeak00]=funcShiftPeak(peak00,rows,columns);
figure(24);
funcPlotAbsImageLineProbe(shiftedPeak00,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%9. Inverse FFT of (0,0) peak
%--------------------------------------------------------------------------
result00Sample=fft2(shiftedPeak00); 
figure(25);
[g1 g2]=funcPlotAbsArrayLineProbe(result10Sample,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%10. Extract (0,1) peak
%--------------------------------------------------------------------------
[peak01,rows01,columns01]=funcExtractPeak01(sampleFFT);
figure(26);
[g1 g2]=funcPlotAbsImageLineProbe(peak01,rows01);
fprintf('Program paused. Press enter to continue.\n');
pause;


%%
%--------------------------------------------------------------------------
%11. Sine Bell filter (0,1) peak
%--------------------------------------------------------------------------
figure(27)
subplot(1,2,1)
gBeforeSineBell=plot(abs(peak01(1+round(rows01/2),:)),'r');
hold on
plot(real(peak01(1+round(rows01/2),:)),'b');
hold on
plot(imag(peak01(1+round(rows01/2),:)),'g');
peak01=funcSineBellFilter(peak01);
subplot(1,2,2)
gAfterSineBell=plot(abs(peak01(1+round(rows01/2),:)),'r');
hold on
plot(real(peak01(1+round(rows01/2),:)),'b');
hold on
plot(imag(peak01(1+round(rows01/2),:)),'g');
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%12. Place (0,1) peak in center of image
%--------------------------------------------------------------------------
[shiftedPeak01]=funcShiftPeak(peak01,rows,columns);
figure(28);
funcPlotAbsImageLineProbe(shiftedPeak01,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%13. Inverse FFT of (0,1) peak
%--------------------------------------------------------------------------
result01Sample=fft2(shiftedPeak01); 
figure(29);
[g1 g2]=funcPlotAbsArrayLineProbe(result01Sample,rows+2);
fprintf('Program paused. Press enter to continue.\n');
pause;


%%
%-------------------------------------------------------------------------------------------------
%Phase:Horizontal
%-------------------------------------------------------------------------------------------------
%flagGraph

%%
%--------------------------------------------------------------------------
%1. Calibrated harmonic image of Peak01
%--------------------------------------------------------------------------
calibratedHarmonicImageHorizontal=result10Sample./result10Ref;

%%
%--------------------------------------------------------------------------
%2. Magnitude image of Peak01
%--------------------------------------------------------------------------
magnitudeHorizontal=abs(calibratedHarmonicImageHorizontal);
figure(30);
[g1 g2]=funcPlotArrayLineProbe(magnitudeHorizontal,rows+2,0,1.5);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%3. Phase image of Peak01
%--------------------------------------------------------------------------
phaseHorizontal=atan2(imag(calibratedHarmonicImageHorizontal),real(calibratedHarmonicImageHorizontal));
figure(31);
[g1 g2]=funcPlotArrayLineProbe(phaseHorizontal,rows+2,-pi,pi);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%-------------------------------------------------------------------------------------------------
%Phase:Vertical
%-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Calibrated harmonic image of Peak01
%--------------------------------------------------------------------------
calibratedHarmonicImageVertical=result01Sample./result01Ref;%A./B is the matrix with elements A(i,j)/B(i,j).

%%
%--------------------------------------------------------------------------
%2. Magnitude image of Peak01
%--------------------------------------------------------------------------
magnitudeVertical=abs(calibratedHarmonicImageVertical);
figure(32);
[~, g2]=funcPlotArrayLineProbe(magnitudeVertical,rows+2,0,1.5);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%--------------------------------------------------------------------------
%3. Phase image of Peak01
%--------------------------------------------------------------------------
phaseVertical=atan2(imag(calibratedHarmonicImageVertical), real(calibratedHarmonicImageVertical));
figure(33);
[g1 g2]=funcPlotArrayLineProbe(phaseVertical,rows+2,-pi/4,pi/4);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%-------------------------------------------------------------------------------------------------
%Phase:Combine Horizontal and Vertical
%-------------------------------------------------------------------------------------------------
fprintf('Dimensions of phaseVertical and phaseHorizontal.\n');
size(phaseVertical)
size(phaseHorizontal)
[rows,columns]=deal(1024,1280);

tempSum=phaseHorizontal+1i*phaseVertical;
numerator=tempSum(1:rows,1:columns);
size(numerator);
numerator=ifft2(numerator);
numerator=circshift(numerator,[0 round(columns/2)]);
numerator=transpose(circshift(transpose(numerator),[0 round(rows/2)]));%numerator=fftshift(numerator);
fprintf('Dimensions of numerator.\n');
size(numerator)
figure(34);
imresize(imadjust(real(imagesc(log10(0.1+abs(numerator))))),256);
colormap(gray)
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%Mathematica code: denominator = Table[ 2 \[Pi] I (k + I l), {k, -1, 1, 2/(rows - 1) // N}, {l, -1, 1, 2/(columns - 1) // N}];
k=-1:2/(rows-1):1;
l=-1:2/(columns-1):1;
denominator=zeros(length(k),length(l));%size:1024x1280
for p=1:length(k); 
    for q=1:length(l);
denominator(p,q)=2*pi*1i*(k(p)+1i*l(q));
    end
end

size(denominator);
[i,j]=find(abs(denominator)==0);
index=[i,j];
if length(index)>0 
    for i=1:length(index);
        r=index(i,1);
        c=index(i,2);
denominator(r,c)=min(2/(rows-1),2/(columns-1));
    end
end
fraction=numerator./denominator;%//N
size(fraction);


phase=abs(ifft2(fraction));
figure(35);
imresize(imadjust(real(imagesc(phase))),256);
colormap(gray)
fprintf('Program paused. Press enter to continue.\n');
pause;

absPhase=abs(phase);
clim=[min(absPhase(:)) max(absPhase(:))];
figure(36);
subplot(1,2,1)
gPhase=imresize(imagesc(phase,clim),500);
title('phase');
colorbar;
colormap(gray);
subplot(1,2,2)
gLinePlot=plot(phase(round(rows/2),:));
fprintf('Program paused. Press enter to continue.\n');
pause;


%%
%-------------------------------------------------------------------------------------------------
%Absorption
%-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Calculate Absorption: Abs = -ln[  | I(0,0)' | ] 
%--------------------------------------------------------------------------
absorption=-log(abs(result00Sample./result00Ref));%matrix division
figure(37);
[g1 g2]=funcPlotArrayLineProbe(absorption,rows+2,-0.2,1.5);
fprintf('Program paused. Press enter to continue.\n');
pause;




%%
%-------------------------------------------------------------------------------------------------
%Diffraction: Horizontal
%-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Calculate Diffraction:  D = -ln[ | I(1,0) | / | I(0,0) | ]    Eq 2 from Wen 2010
%--------------------------------------------------------------------------
diffractionHorizontal=-log(abs(result10Sample./result10Ref)./abs(result00Sample./result00Ref));%log or log10?
figure(38);
[g1 g2]=funcPlotArrayLineProbe(diffractionHorizontal,rows+2,0,2);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%-------------------------------------------------------------------------------------------------
%Diffraction: Vertical
%-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Calculate Diffraction:  D = -ln[ | I(1,0) | / | I(0,0) | ]    Eq 2 from Wen 2010
%--------------------------------------------------------------------------
diffractionVertical=-log(abs(result01Sample./result01Ref)./abs(result00Sample./result00Ref));
figure(39);
[g1 g2]=funcPlotArrayLineProbe(diffractionVertical,rows+2,0,6);
fprintf('Program paused. Press enter to continue.\n');
pause;

%%
%-------------------------------------------------------------------------------------------------
%Diffraction: Combine Horizontal and Vertical
%-------------------------------------------------------------------------------------------------
diffraction=abs(diffractionVertical+1i*diffractionHorizontal);
[min(diffraction(:)),max(diffraction(:))];
figure(40);
imresize(imadjust(imagesc(real(diffraction))),256); 
colormap(gray)
fprintf('Program paused. Press enter to continue.\n');
pause;

clim=[min(diffraction(:)),max(diffraction(:))];
figure(41);
gPhase=imresize(imagesc(diffraction,clim),500);
colorbar;
colormap(gray)

fprintf('Do you want to delete all the figures? If so, please press enter.\n');
pause;
%delete figures
all_figs = findall(0, 'type', 'figure');
delete(all_figs);
break
%%
%-------------------------------------------------------------------------------------------------
%Save absorption, phase(3), and dark field(3)
%-------------------------------------------------------------------------------------------------
%haven't done


%%
%Process all Sample files
%Exploration: Use larger zero-pad, IFFT, crop: Will this remove truncation wiggles?
%-------------------------------------------------------------------------------------------------
%old code
%-------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------
%12. Place (0,1) peak in center of image
%--------------------------------------------------------------------------
%shiftedPeak01=funcShiftPeak(peak01,rows,columns);
%funcPlotAbsImageLineProbe(shiftedPeak01,rows+2);



%--------------------------------------------------------------------------
%13. Inverse FFT of (0,1) peak
%--------------------------------------------------------------------------
%result01Sample=ifft(shiftedPeak01);
%funcPlotAbsArrayLineProbe(result01Sample,rows+2);


%-------------------------------------------------------------------------------------------------
%new code
%-------------------------------------------------------------------------------------------------
%shiftedPeak01=funcShiftPeakZeroPad(peak01,rows,columns);
%[rows columns]
%size(peak01)
%size(shiftedPeak01)
%result01Sample=funcIFFTandCrop(shiftedPeak01,rows,columns);
%size(result01Sample)
%flagGraph
%funcPlotAbsArrayLineProbe(result01Sample,floor(1.1*rows)+2);



%delete figures
all_figs = findall(0, 'type', 'figure');
delete(all_figs);


%figures to keep
figs2keep = [4, 7];
% Uncomment the following to 
% include ALL windows, including those with hidden handles (e.g. GUIs)
% all_figs = findall(0, 'type', 'figure');
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));
