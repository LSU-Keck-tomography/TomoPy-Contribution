%% Foraminifera
%Melete04 kernel (default kernel for this notebook)
%transmission:a0/a0 (ref)
%absorption:-Log[a0/a0(ref)]
%differential phase: phase-phase(ref)
%darkfield: a1/a1(ref) / a0/a0(ref)

%% Note
%  Please change homeDirectory=cd('/melete-nas01/lsu/scratch/jumao/matlab_stepped-grating')
%                to your home directory, where the Foraminifera data set is

%% Initialization
clear all; clc;
fprintf('Home directory is: \n');
homeDirectory=cd('/melete-nas01/lsu/scratch/jumao/Matlab_stepped-grating');
fprintf('%s \n\n',homeDirectory);
fprintf('Program paused. Press enter to continue. \n');
fprintf('\n');
pause;

%parpool

%%
%-------------------------------------------------------------------------------------------------
%Step 1: Initialize Paths
%-------------------------------------------------------------------------------------------------
pathReferences=[homeDirectory '/projections/'];
pathProjections=pathReferences;
pathSinoAbs=[homeDirectory '/sino_abs/'];
pathSinoDIhorizontal=[homeDirectory '/sino_DI_H/'];
pathSinoDPChorizontal=[homeDirectory '/sino_DPC_H/'];
pathSinoDIvertical=[homeDirectory '/sino_DI_V/'];
pathSinoDPCvertical=[homeDirectory '/sino_DPC_V/'];
pathFigures=[homeDirectory '/figures/'];
pathSlices=[homeDirectory '/slices/'];
nameSample='proj_';


%%
%-------------------------------------------------------------------------------------------------
%Step 2: Create Directories
%-------------------------------------------------------------------------------------------------
fprintf('Creating directories: \n projections \n sino_abs \n sino_DI_H \n sino_DPC_H \n sino_DI_V \n sino_DPC_V \n figures \n slices \n \n');
mkdir(pathReferences)
%mkdir(pathProjections)
mkdir(pathSinoAbs)
mkdir(pathSinoDIhorizontal)
mkdir(pathSinoDPChorizontal)
mkdir(pathSinoDIvertical)
mkdir(pathSinoDPCvertical)
mkdir(pathFigures)
mkdir(pathSlices)
fprintf('Program paused. Press enter to continue. \n');
fprintf('\n');
pause;
%%
%-------------------------------------------------------------------------------------------------
%Step 3: Create list of angles
%-------------------------------------------------------------------------------------------------
angleStart=0;
angleLastStop=180;
angleStep=0.3;

angleListFull=(angleStart:angleStep:angleLastStop);
fprintf('Number of angles.\n');
fprintf('%d \n\n', length(angleListFull));

global NZ
numberAngles=(angleLastStop-angleStart)/angleStep+1;
NZ=numberAngles;
fprintf('Program paused. Press enter to continue.\n\n');
pause

%%
%-------------------------------------------------------------------------------------------------
%Step 4. Make all sinograms (Abs, DI, DPC)
%-------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------
%Step 4.0 Define functions for air normalization, streak removal (FFT-high
%pass), sinogram shift, and phase-unwrap (Haas cost function method)
%-------------------------------------------------------------------------------------------------

%funcUseAirRegionForOffset
%funcStreakRemovalFFT
%funcShiftSinogram
%funcUseAirRegionForOffsetABSandDPC
%funchaasPhaseSinogram

%%
%-------------------------------------------------------------------------------------------------
%Step 4.1 Read all projections and make Absorption allSinograms 
%-------------------------------------------------------------------------------------------------
filenamesProjection=dir([pathProjections nameSample '*.h5']);
fprintf('The filename of the first projection file is:\n');
fprintf('%s \n\n',filenamesProjection(1).name);
fprintf('Number of projections.\n\n');
fprintf('%d \n\n',length(filenamesProjection));

firstDataAbsorption=h5read([pathProjections filenamesProjection(1).name],'/absorption');
[numberSlices numberColumns]=size(transpose(firstDataAbsorption));
[NY NX]=deal(numberSlices,numberColumns);
numberRays=numberColumns;

fprintf('Dimensions of all projections.\n');
fprintf('%d \n\n', [NX NY NZ]);
fprintf('Program paused. Press enter to continue.\n\n');
pause

%Show memory in use
allSinograms=zeros(NZ,NX,NY);
fprintf('Dimensions of allSinograms.\n');
fprintf('%f \n\n', size(allSinograms))
fprintf('Memory used:\n');
unix('vm_stat')
fprintf('Program paused. Press enter to continue.\n\n');
pause

%Initial time before loop
t0=cputime;

fprintf('Loop evaluating.\n');
for indexProjection=1:NZ
filename=filenamesProjection(indexProjection);
dataAbsorption=h5read([pathProjections filename.name],'/absorption');
allSinograms(indexProjection,:,:)=dataAbsorption;
end
allSinograms=permute(allSinograms, [2, 1, 3]);

fprintf('Writing absorption binary file.\n\n');
filename = [pathSinoAbs 'sino_all_abs_Matlab.bin'];
fid=fopen(filename,'wb');
fwrite(fid,allSinograms,'real*4');

%Final time after loop; total time loop ran (measure of performance)
t1=cputime;
dateDifference = (t1 - t0)/60;
fprintf('Approximate loop CPU time: \n');
fprintf('%f \n\n', dateDifference)

%Show sinogram with line probe
index=165;%Any number from 1 to 520
oneSinogram=allSinograms(:,:,index);
[min(oneSinogram(:)) max(oneSinogram(:))];
figure(1);
subplot(1,2,1)
gSino=imagesc(oneSinogram);
colorbar;
colormap(gray)

lineProbe=allSinograms([30 120],:,index);
lineProbe2D=reshape(lineProbe,2,NZ);
subplot(1,2,2)
for i=1:2
gLineProbe=plot(lineProbe2D(i,:));
hold on
end
fprintf('Program paused. Press enter to continue.\n\n');
pause
%%
%-------------------------------------------------------------------------------------------------
%Step 4.2 Read all projections and make DI(Vertical) allSinograms 
%-------------------------------------------------------------------------------------------------
filenamesProjection=dir([pathProjections nameSample '*.h5']);
fprintf('The filename of the first projection file is:\n');
fprintf('%s \n\n', filenamesProjection(1).name);
fprintf('Number of projections.\n');
fprintf('%d \n\n', length(filenamesProjection));

%Show memory in use
allSinograms=zeros(NZ,NX,NY);
fprintf('Dimensions of allSinograms.\n');
fprintf('%d \n\n', size(allSinograms));
fprintf('Memory used:\n');
unix('vm_stat')
fprintf('Program paused. Press enter to continue.\n\n');
pause

%Initial time before loop
t0 = cputime;

fprintf('Loop evaluating.\n');
for indexProjection=1:NZ
filename=filenamesProjection(indexProjection);
dataAbsorption=h5read([pathProjections filename.name],'/diffraction');
allSinograms(indexProjection,:,:)=dataAbsorption;
end
allSinograms=permute(allSinograms, [2, 1, 3]);

fprintf('Writing dark field binary file.\n\n');
filename = [pathSinoDIvertical 'sino_all_DI_vertical_Matlab.bin'];
fid=fopen(filename,'wb');
fwrite(fid,allSinograms,'real*4');

%Final time after loop; total time loop ran (measure of performance)
t1 = cputime;
dateDifference = (t1 - t0)/60;
fprintf('Approximate loop CPU time: \n');
fprintf('%f \n\n', dateDifference)

%Show sinogram with line probe
index=165;%Any number from 1 to 520
oneSinogram=allSinograms(:,:,index);
[min(oneSinogram(:)) max(oneSinogram(:))];
figure(2);
subplot(1,2,1)
gSino=imagesc(oneSinogram);
colorbar;
colormap(gray)

lineProbe=allSinograms([30 120],:,index);
lineProbe2D=reshape(lineProbe,2,NZ);
subplot(1,2,2)
for i=1:2
gLineProbe=plot(lineProbe2D(i,:));
hold on
end
fprintf('Program paused. Press enter to continue.\n\n');
pause
%%
%-------------------------------------------------------------------------------------------------
%Step 4.3 Read all projections and make DPC(Vertical) allSinograms 
%-------------------------------------------------------------------------------------------------
filenamesProjection=dir([pathProjections nameSample '*.h5']);
fprintf('The filename of the first projection file is:\n');
fprintf('%s \n\n', filenamesProjection(1).name);
fprintf('Number of projections.\n');
fprintf('%d \n\n', length(filenamesProjection))

%Show memory in use
allSinograms=zeros(NZ,NX,NY);
fprintf('Dimensions of allSinograms.\n');
fprintf('%d \n\n', size(allSinograms));
fprintf('Memory used:\n');
unix('vm_stat')
fprintf('Program paused. Press enter to continue.\n\n');
pause

%Initial time before loop
t0 = cputime;

fprintf('Loop evaluating.\n');
for indexProjection=1:NZ
filename=filenamesProjection(indexProjection);
dataAbsorption=h5read([pathProjections filename.name],'/differentialPhase');
allSinograms(indexProjection,:,:)=dataAbsorption;
end
allSinograms=permute(allSinograms, [2, 1, 3]);

fprintf('Writing DPC binary file.\n\n');
filename = [pathSinoDPCvertical 'sino_all_DPC_vertical_Matlab.bin'];
fid=fopen(filename,'wb');
fwrite(fid,allSinograms,'real*4');

%Final time after loop; total time loop ran (measure of performance)
t1 = cputime;
dateDifference = (t1 - t0)/60;
fprintf('Approximate loop CPU time: \n');
fprintf('%f \n\n', dateDifference)

%Show sinogram with line probe
index=165;%Any number from 1 to 520
oneSinogram=allSinograms(:,:,index);
[min(oneSinogram(:)) max(oneSinogram(:))];
figure(3);
subplot(1,2,1)
gSino=imagesc(oneSinogram);
colorbar;
colormap(gray)

lineProbe=allSinograms([30 120],:,index);
lineProbe2D=reshape(lineProbe,2,NZ);
subplot(1,2,2)
for i=1:2
gLineProbe=plot(lineProbe2D(i,:));
hold on
end
fprintf('Program paused. Press enter to continue.\n\n');
pause

%%
%-------------------------------------------------------------------------------------------------
%Step 4.5 Remove streaks, apply air offset, shift the sinograms, and save
%as HDF5 files
%-------------------------------------------------------------------------------------------------
fprintf('Remove streaks, apply air offset, shift sinograms, and save as HDF5 files. \nPress enter to continue.\n');
pause;
clear allSinograms
bestShift=-12;

%% Absorption HDF5
filename=[pathSinoAbs 'sino_all_abs_Matlab.bin'];
fid=fopen(filename,'rb');
allSinograms=fread(fid,'real*4');

fprintf('Dimensions of allSinograms.\n');
fprintf('%d \n\n', size(allSinograms));
allSinograms=reshape(allSinograms,[NX NZ NY]);
fprintf('Dimensions of allSinograms after reshape.\n');
fprintf('%d \n\n', size(allSinograms));
fprintf('Program paused. Press enter to continue.\n\n');
pause

fprintf('Writing absorption sinograms.\n');
fprintf('Press enter to continue\n\n');
for indexSlice=1:numberSlices
sinogram=allSinograms(:,:,indexSlice);
sinogram=funcUseAirRegionForOffset(sinogram,10);
sinogram=funcShiftSinogram(sinogram,bestShift,numberAngles,numberColumns);
sinogram=transpose(sinogram);
filenameHDF5=[pathSinoAbs 'sino_abs_' num2str(indexSlice,'%0.4d') '.h5'];
[NNX NNY]=size(sinogram);
h5create(fullfile(filenameHDF5),'/sinogram',[NNX NNY]);
h5write(fullfile(filenameHDF5), '/sinogram',sinogram);
end

%% Dark Field HDF5
filename=[pathSinoDIvertical 'sino_all_DI_vertical_Matlab.bin'];
fid=fopen(filename,'rb');
allSinograms=fread(fid,'real*4');
fprintf('Dimensions of allSinograms.\n');
fprintf('%d \n\n', size(allSinograms));
allSinograms=reshape(allSinograms,[NX NZ NY]);
fprintf('Dimensions of allSinograms after reshape.\n');
fprintf('%d \n\n', size(allSinograms));
fprintf('Program paused. Press enter to continue.\n\n');

fprintf('Writing dark field sinograms.\n');
fprintf('Press enter to continue\n\n');
for indexSlice=1:numberSlices
sinogram=allSinograms(:,:,indexSlice);
sinogram=1-sinogram;
sinogram=funcUseAirRegionForOffset(sinogram,10);
%sinogram=funcStreakRemovalFFT(sinogram,0.005);
sinogram=funcShiftSinogram(sinogram,bestShift,numberAngles,numberColumns);
sinogram=transpose(sinogram);
filenameHDF5=[pathSinoDIvertical 'sino_DI_V_' num2str(indexSlice,'%0.4d') '.h5'];
[NNX NNY]=size(sinogram);
h5create(fullfile(filenameHDF5),'/sinogram',[NNX NNY]);
h5write(fullfile(filenameHDF5), '/sinogram',sinogram);
end
%% DPC HDF5
fprintf('Writing DPC sinograms.\n');
fprintf('Press enter to continue\n\n');
filename=[pathSinoDPCvertical 'sino_all_DPC_vertical_Matlab.bin'];


for indexSlice=1:numberSlices
fid=fopen(filename,'rb');
fseek(fid, 4*NX*NZ*(indexSlice-1), 'bof');
sinogram=fread(fid, NX*NZ,'real*4');
sinogram=reshape(sinogram,[NX, NZ]);
sinogram=funcUseAirRegionForOffset(sinogram,10);
%sinogram=funcStreakRemovalFFT(sinogram,0.001);
sinogram=funcShiftSinogram(sinogram,bestShift,numberAngles,numberColumns);
%filenameABS=[pathSinoAbs 'sino_abs_' num2str(indexSlice,'%0.4d') '.h5'];
%sinogramABS=h5read(filenameABS,'/sinogram');
%sinogram=funchaasPhaseSinogram(sinogram,sinogramABS,NX,NZ);
sinogram=transpose(sinogram);
filenameHDF5=[pathSinoDPCvertical 'sino_DPC_V_' num2str(indexSlice,'%0.4d') '.h5'];
[NNX NNY]=size(sinogram);
h5create(filenameHDF5,'/sinogram',[NNX NNY]);
h5write(filenameHDF5, '/sinogram',sinogram);
end

%%
%-------------------------------------------------------------------------------------------------
%Step 4.6 Find best center of rotation using funcShiftSinogram and reconstruction of a slice: 
%bestShifted=-12 
%-------------------------------------------------------------------------------------------------
filename =[pathSinoAbs 'sino_all_abs_Matlab.bin'];
fid = fopen(filename,'rb');
allSinograms = fread(fid,'real*4');
fprintf('Dimensions of allSinograms.\n');
fprintf('%d \n\n', size(allSinograms));

allSinograms=reshape(allSinograms,[NX NZ NY]);
fprintf('Dimensions of allSinograms after reshape.\n');
fprintf('%d \n\n', size(allSinograms));
[numberColumns, numberAngles, numberSlices] = size(allSinograms);

index=400;%Any number from 1 to 520
oneSinogram=allSinograms(:,:,index);
fprintf('Dimensions of oneSinogram.\n');
fprintf('%d \n\n', size(oneSinogram));

%Sinogram
figure(4);
subplot(3,1,1)
gSino=imagesc(oneSinogram);
colorbar;
colormap(gray)

%Line Probe
subplot(3,1,2);
lineProbe=allSinograms(30:round(numberColumns/2),:,index);
lineProbe=reshape(lineProbe,round(numberColumns/2)-30+1,NZ);
for i=1:round(numberColumns/2)+1-30
gLineProbe=plot(lineProbe(i,:));
hold on
end

%Sinogram
subplot(3,1,3);
gThumb=imresize(imadjust(imagesc(real(oneSinogram))),700);
fprintf('Program paused. Press enter to continue.\n\n');

%Image Data
sinogram=oneSinogram;
sinogram=funcUseAirRegionForOffset(sinogram,10);
%sinogram=funcStreakRemovalFFT(sinogram,0.005);
sinogram=sinogram(:,(1:2:numberAngles));
[newNumberRays,newNumberAngles]=size(sinogram);

R=real(funcShiftSinogram(sinogram, bestShift,newNumberAngles, newNumberRays));
fprintf('Dimensions of shifted sinogram\n')
fprintf('%d \n\n', size(R));

%Slice
figure(5);
irR=iradon(R, 0:0.6:180,'linear','Cosine',1.0,numberRays);
imageR=imagesc(irR);
imadjust(imresize(imageR,600));
colorbar;
colormap(gray);
%%
%-------------------------------------------------------------------------------------------------
%Step 4.7 Clear variables
%-------------------------------------------------------------------------------------------------
fprintf('Clear allSinograms.\n');
clear allSinograms;

%%
%-------------------------------------------------------------------------------------------------
%Step 5. Reconstruct the sinograms into slices (Prefer Matlab)
%-------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------
%Step 5.1 Read one sinograms and reconstruct into a slice. Test different
%kernels(poor performance for most kernels, TS 36054 Oct8,2012)
%-------------------------------------------------------------------------------------------------
%Use iradon to build slices
fprintf('Read one sinograms and reconstruct into a slice.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');

clear numberSlices
clear numberAngles
filenamesSinogram=[pathSinoAbs 'sino_abs_0100.h5'];
sinogram=h5read(filenamesSinogram, '/sinogram');
sinogram=transpose(sinogram);
[numberAngles numberRays]=size(sinogram);

%Various filters
sliceRamp=iradon(sinogram,angleListFull,'linear','Ram-Lak',1.0,numberRays);
sliceRampCosine=iradon(sinogram,angleListFull,'linear','Cosine',1.0,numberRays);
sliceHann=iradon(sinogram,angleListFull,'linear','Hann',1.0,numberRays);
sliceButterHamming=iradon(sinogram,angleListFull,'linear','Hamming',1.0,numberRays);
sliceSheppLogan=iradon(sinogram,angleListFull,'linear','Shepp-Logan',1.0,numberRays);
sliceNone=iradon(sinogram,angleListFull,'linear','None',1.0,numberRays);

%Slice filtered images for comparison
fprintf('Show images of each slice.\n\n');
figure(6);
subplot(3,2,1);
imadjust(imresize(imagesc(sliceRamp),700));
subplot(3,2,2);
imadjust(imresize(imagesc(sliceRampCosine),700));
subplot(3,2,3);
imadjust(imresize(imagesc(sliceHann),700));
subplot(3,2,4);
imadjust(imresize(imagesc(sliceButterHamming),700));
subplot(3,2,5);
imadjust(imresize(imagesc(sliceSheppLogan),700));
subplot(3,2,6);
imadjust(imresize(imagesc(sliceNone),700));

colormap(gray)
%%
%-------------------------------------------------------------------------------------------------
%Step 5.2 Absorption:  Read all sinograms and reconstruct into slices
%-------------------------------------------------------------------------------------------------
clear numberSlices
fprintf('Absorption:  Read all sinograms and reconstruct into slices.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filenamesSinograms=dir([pathSinoAbs 'sino_abs_' '*.h5']);
numberSlices=length(filenamesSinograms);
filenamesSinogramsOne=[pathSinoAbs filenamesSinograms(1).name];
sinogram=h5read(filenamesSinogramsOne,'/sinogram');
sinogram=transpose(sinogram);
[numberColumns,numberRows]=size(sinogram);

for indexSlice=1:numberSlices
    filenamesSinogramsIndex=[pathSinoAbs filenamesSinograms(indexSlice).name];
    sinogram=h5read(filenamesSinogramsIndex,'/sinogram');
    sinogram=transpose(sinogram);
    slice=iradon(sinogram,angleListFull,'linear','Cosine',1.0,numberRays);
    filenameHDF5=[pathSlices 'AbsSlice_' num2str(indexSlice,'%0.4d') '.h5'];
    [NNX NNY]=size(slice);
    h5create(filenameHDF5, '/absslice', [NNX NNY]);
    h5write(filenameHDF5, '/absslice', slice);
end

%%
%-------------------------------------------------------------------------------------------------
%Step 5.3 Dark Field:  Read all sinograms and reconstruct into slices
%-------------------------------------------------------------------------------------------------
clear numberSlices
fprintf('Dark Field:  Read all sinograms and reconstruct into slices.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filenamesSinograms=dir([pathSinoDIvertical 'sino_DI_V_' '*.h5']);
numberSlices=length(filenamesSinograms);
filenamesSinogramsOne=[pathSinoDIvertical filenamesSinograms(1).name];
sinogram=h5read(filenamesSinogramsOne,'/sinogram');
sinogram=transpose(sinogram);
[numberColumns,numberRows]=size(sinogram);

for indexSlice=1:numberSlices
    filenamesSinogramsIndex=[pathSinoDIvertical filenamesSinograms(indexSlice).name];
    sinogram=h5read(filenamesSinogramsIndex,'/sinogram');
    sinogram=transpose(sinogram);
    slice=iradon(sinogram,angleListFull,'linear','Cosine',1.0,numberRays);
    filenameHDF5=[pathSlices 'DISlice_' num2str(indexSlice,'%0.4d') '.h5'];
    [NNX NNY]=size(slice);
    h5create(filenameHDF5,'/dislice',[NNX NNY]);
    h5write(filenameHDF5, '/dislice', slice);
end

%%
%-------------------------------------------------------------------------------------------------
%Step 5.4 DPC:  Read all sinograms and reconstruct into slices
%-------------------------------------------------------------------------------------------------
clear numberSlices
fprintf('DPC:  Read all sinograms and reconstruct into slices.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filenamesSinograms=dir([pathSinoDPCvertical 'sino_DPC_V_' '*.h5']);
numberSlices=length(filenamesSinograms);
filenamesSinogramsOne=[pathSinoDPCvertical filenamesSinograms(1).name];
sinogram=h5read(filenamesSinogramsOne,'/sinogram');
sinogram=transpose(sinogram);
[numberColumns,numberRows]=size(sinogram);

for indexSlice=1:numberSlices
    filenamesSinogramsIndex=[pathSinoDPCvertical filenamesSinograms(indexSlice).name];
    sinogram=h5read(filenamesSinogramsIndex,'/sinogram');
    sinogram=transpose(sinogram);
    slice=iradon(sinogram,angleListFull,'linear','Cosine',1.0,numberRays);
    filenameHDF5=[pathSlices 'DPCSlice_' num2str(indexSlice,'%0.4d') '.h5'];
    [NNX NNY]=size(slice);
    h5create(filenameHDF5,'/dpcslice',[NNX NNY]);
    h5write(filenameHDF5, '/dpcslice', slice);
end
%%
%-------------------------------------------------------------------------------------------------
%Step 5.5 Look at one slice in a sinogram (Abs)
%-------------------------------------------------------------------------------------------------
fprintf('Look at one slice in a sinogram (Abs).\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filename=[pathSinoAbs 'sino_all_abs_Matlab.bin'];
indexSlice=502; %Any number from 1 to 520
fid=fopen(filename,'rb');
fseek(fid, 4*NX*NZ*(indexSlice-1), 'bof');
sinogram=fread(fid, NX*NZ,'real*4');
sinogram=reshape(sinogram,[NX, NZ]);
fprintf('Dimension of sinogram.\n\n');
fprintf('%d \n\n', size(sinogram));
figure(7);
imresize(imadjust(imagesc(real(sinogram))),300);

%%
%-------------------------------------------------------------------------------------------------
%Step 6. Store all slices in one binary file
%-------------------------------------------------------------------------------------------------
%Scale factors
fprintf('Store all slices in one binary file.\n\n');
scaleFactorAbs=3000000;
scaleFactorDI=1000000;
scaleFactorDPC=10000;
nameSample;

%%
%-------------------------------------------------------------------------------------------------
%Step 6.1 Absorption
%-------------------------------------------------------------------------------------------------
fprintf('Absorption.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filenamesSlices=dir([pathSlices 'AbsSlice_' '*.h5']);
numberSlices=length(filenamesSlices);
filenamesSlicesOne=[pathSlices filenamesSlices(1).name];
slice=h5read(filenamesSlicesOne,'/absslice');
[numberRays,numberColumns]=size(slice);
volume=zeros([numberColumns,numberColumns,numberSlices]);

%Scale factor options for filenameRAW
if scaleFactorAbs == 5000
    strScale = '_5k';
elseif scaleFactorAbs == 100000
    strScale = '_100k';
elseif scaleFactorAbs == 1000000
    strScale = '_1M';
else
    strScale = '';
end
filenameRAW=[pathSlices nameSample 'abs_' num2str(numberColumns) '_' num2str(numberColumns) '_' num2str(numberSlices) strScale '_int16_le.raw'];

%Initial time before loop
t0 = cputime;

for indexSlice=1:numberSlices;
    filenamesSlicesIndex=[pathSlices filenamesSlices(indexSlice).name];
    sliceReal=h5read(filenamesSlicesIndex,'/absslice');
    sliceReal=transpose(sliceReal);
    sliceInteger=round(sliceReal*scaleFactorAbs);
    large=find(sliceInteger > (2^15 - 1));
    small=find(sliceInteger < (-2^15));
    sliceInteger(large) = 2^15 - 1;
    sliceInteger(small) = -2^15;
    volume(:,:,indexSlice)=sliceInteger;
end

fid=fopen(filenameRAW,'wb');
fwrite(fid,volume,'int16');

%Final time after loop; total time loop ran (measure of performance)
t1 = cputime;
dateDifference = (t1 - t0)/60;
fprintf('Approximate loop CPU time: \n');
fprintf('%f \n\n', dateDifference);

%%
%-------------------------------------------------------------------------------------------------
%Step 6.2 Dark Field
%-------------------------------------------------------------------------------------------------
fprintf('Dark Field.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filenamesSlices=dir([pathSlices 'DISlice_' '*.h5']);
numberSlices=length(filenamesSlices);
filenamesSlicesOne=[pathSlices filenamesSlices(1).name];
slice=h5read(filenamesSlicesOne,'/dislice');
[numberRays,numberColumns]=size(slice);
volume=zeros([numberColumns,numberColumns,numberSlices]);

%Scale factor options for filenameRAW
if scaleFactorDI == 5000
    strScale = '_5k';
elseif scaleFactorDI == 100000
    strScale = '_100k';
elseif scaleFactorDI == 1000000
    strScale = '_1M';
else
    strScale = '_';
end
filenameRAW=[pathSlices nameSample 'di_' num2str(numberColumns) '_' num2str(numberColumns) '_' num2str(numberSlices) strScale '_int16_le.raw'];

%Initial time before loop
t0=cputime;

for indexSlice=1:numberSlices;
    filenamesSlicesIndex=[pathSlices filenamesSlices(indexSlice).name];
    sliceReal=h5read(filenamesSlicesIndex,'/dislice');
    sliceReal=transpose(sliceReal);
    sliceInteger=round(sliceReal*scaleFactorAbs);
    large=find(sliceInteger > (2^15 - 1));
    small=find(sliceInteger < (-2^15));
    sliceInteger(large) = 2^15 - 1;
    sliceInteger(small) = -2^15;
    volume(:,:,indexSlice)=sliceInteger;
end

fid=fopen(filenameRAW,'wb');
fwrite(fid,volume,'int16');

%Final time after loop; total time loop ran (measure of performance)
t1=cputime;
dateDifference = (t1 - t0)/60;
fprintf('Approximate loop CPU time: \n');
fprintf('%f \n\n', dateDifference);

%%
%-------------------------------------------------------------------------------------------------
%Step 6.3 DPC
%-------------------------------------------------------------------------------------------------
fprintf('Differential Phase Contrast.\n\n');
fprintf('Program paused. Press enter to continue.\n\n');
filenamesSlices=dir([pathSlices 'DPCSlice_' '*.h5']);
numberSlices=length(filenamesSlices);
filenamesSlicesOne=[pathSlices filenamesSlices(1).name];
slice=h5read(filenamesSlicesOne,'/dpcslice');
[numberRays,numberColumns]=size(slice);
volume=zeros([numberColumns,numberColumns,numberSlices]);

%Scale factor options for filenameRAW
if scaleFactorDPC == 5000
    strScale = '_5k';
elseif scaleFactorDPC == 100000
    strScale = '_100k';
elseif scaleFactorDPC == 1000000
    strScale = '_1M';
else
    strScale = '_';
end
filenameRAW=[pathSlices nameSample 'dpc_' num2str(numberColumns) '_' num2str(numberColumns) '_' num2str(numberSlices) strScale '_int16_le.raw'];

%Initial time before loop
t0=cputime;

for indexSlice=1:numberSlices;
    filenamesSlicesIndex=[pathSlices filenamesSlices(indexSlice).name];
    sliceReal=h5read(filenamesSlicesIndex,'/dpcslice');
    sliceReal=transpose(sliceReal);
    sliceInteger=round(sliceReal*scaleFactorAbs);
    large=find(sliceInteger > (2^15 - 1));
    small=find(sliceInteger < (-2^15));
    sliceInteger(large) = 2^15 - 1;
    sliceInteger(small) = -2^15;
    volume(:,:,indexSlice)=sliceInteger;
end

fid=fopen(filenameRAW,'wb');
fwrite(fid,volume,'int16');

%Final time after loop; total time loop ran (measure of performance)
t1=cputime;
dateDifference = (t1 - t0)/60;
fprintf('Approximate loop CPU time: \n');
fprintf('%f \n\n', dateDifference);

%% delete figures
fprintf('Do you want to delete all the figures? If so, please press enter.\n');
pause
all_figs = findall(0, 'type', 'figure');
delete(all_figs);