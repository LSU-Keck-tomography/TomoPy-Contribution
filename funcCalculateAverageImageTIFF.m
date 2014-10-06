function [dataAverage]=funcCalculateAverageImageTIFF(filenames)
for k=1:length(filenames)         
filenamesData{k}=imread(filenames(k).name);
end
%[rows,columns]=size(filenamesData{k});

global rows columns
[rows,columns]=size(filenamesData{1});


dataAverage=zeros([rows,columns]);
for index=1:length(filenames)
dataAverage=double(dataAverage)+double(imread(filenames(index).name));
end
dataAverage=round(dataAverage/length(filenames));
%size(dataAverage);
