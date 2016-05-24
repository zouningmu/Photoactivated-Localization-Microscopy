function iQ_mvtoimg(folderName)
%% SAVE MV AS SEPARATED IMAGE FILES
% INPUT: folderName, have to be a string
tic
cd(folderName); D=dir;
num_file=length(D);
mkdir('images');

total_frame=0;
for j=3:num_file
    fileName=D(j).name;
    tiffInfo = imfinfo(fileName);  %# Get the TIFF file information
    no_frame = numel(tiffInfo);    %# Get the number of images in the file
    total_frame=no_frame+total_frame;
end

name=fullnum(total_frame);        %# Generates filename index

no_frame_acc=0; ind_frame=cell(num_file,1); ind_frame{1,1}=0;
for j=3:num_file
    fileName=D(j).name;
    tiffInfo = imfinfo(fileName);  %# Get the TIFF file information
    no_frame = numel(tiffInfo);    %# Get the number of images in the file
    no_frame_acc=no_frame+no_frame_acc;
    ind_frame{j-1,1}=no_frame_acc;
    heading =cell2mat(regexp(fileName,'(\w+).tif','tokens','once')) ;
   for iFrame = 1:no_frame
     Movie = double(imread(fileName,'Index',iFrame,'Info',tiffInfo));
     Movie=uint16(Movie);
     newname=strcat('no',name{iFrame+ind_frame{j-2,1}},'_',heading);
     cd images
     imwrite(Movie,newname,'tiff');
     cd ..
   end
end
toc
%% nested function fullnum
function b=fullnum(n)
% this function generates series number strings with 0 heading according to the input n.
% for example, if you have 4000 files (n=4000) to save, this function will
% generate a column of strings from 0001 to 4000. The first file with be 0001
% the 100th file will be 0100, and the 1500th file will be 1500.

if n~=floor(n)||n<=0||n>999999
  fprintf('Error:The number is not allowed float-point,less than 0 or larger than 9999!\n');
  return;
end
b={};
for i=1:n
  if i<10^floor(log10(n))
    b{i}=strcat(char(ones(1,floor(log10(n))-floor(log10(i)))*48),int2str(i));
  else
    b{i}=int2str(i);
  end
end
end
end