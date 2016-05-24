function [im_bk,peak,intensity,bk]=iQ_pkfnd(im,threshold,sz)
%% THIS FUNCTION IS FOR FINDING EACH BRIGHT OBJECT IN THE FIELD OF VIEW
% INPUTS:
%       im: image files
%       threshold: define the minimum values for a peak, discard all local
%                  maiximum peaks with intensities lower than this value,
%                  suggested value for the threshold is 4. 
%       sz: intensity window size, suggested value for sz is 11.          
% OUTPUTS:
% im_bk: background removed image
% peak: coordinate list of local maximum spots
% intensity: peak intensity (intensity sum over the mask size)
% bk: [[im_bk_mean,im_bk_std,th] which is [ mean value of im_bk, standard dev of im_bk, threshold value(im_bk_mean+threshold*im_bk_std)]
%%
if nargin < 3, sz = 11; end
if nargin < 2, threshold = 4; end

[nr,nc]=size(im);                                   % find out the size of image
bk_fit = imopen(im,strel('disk',15));               % fit the background
im_bk=im-bk_fit;                                    % remove the background from original image

mask=ones(sz);                                      % define the mask for one PSF spot   
mask(ceil(sz/2),ceil(sz/2))=0;   

bw=im_bk>imdilate(im_bk,mask);                      % find out the local maximum spot with specfied box area

[r,c]=find(bw>0);  rc=[r,c];                        % get the coordinates of local maximum spots
                                           
im_bk_mean=mean(mean(im_bk));                       % find the average value of background removed image
im_bk_std=mean(std(im_bk));                         % find the standard deviation of background removed image
th=im_bk_mean+threshold*im_bk_std;                  % set up the threshold for peak determination
bk=[im_bk_mean,im_bk_std,th];                       % generate the bk info: [im_bk_mean,im_bk_std,th]=[average of image,standard deviation of image,threshold for peak finding

% if nothing is found above the threshold
if isempty(rc)
    display('no local maximum');
    peak=[]; intensity=[];
    return;
else
    rc_th=[];
    r1=1;                                               % obtain the coordinates which shows intensity larger than the threshold
    [rc_row,~]=size(rc);
for q=1:rc_row;
    if im_bk(rc(q,1),rc(q,2))>th;
        rc_th(r1,:)=rc(q,:);
        r1=r1+1;
    end
end
if isempty(rc_th)
    display('nothing above threshold');
    peak=[]; intensity=[];
    return;
else
peak(:,1)=rc_th(:,2);                               % generate the peak output: [column1,column2]= [x,y]
peak(:,2)=rc_th(:,1);

% take spots whose mask is not connected to the boundary, discard the rest spots
ind_boundary= peak(:,1)>ceil(sz/2) & peak(:,1)<nc-ceil(sz/2) & peak(:,2)>ceil(sz/2) & peak(:,2)<nr-ceil(sz/2);
peak=peak(ind_boundary,:); 

ind_th=sub2ind(size(im_bk),peak(:,2),peak(:,1));    % convert new coordinates to index
im_bk_sums=imfilter(im_bk,ones(sz));                % calculate the intergrated intensity over specified area
im_bk_int=im_bk_sums(ind_th);
intensity=im_bk_int;                
end
end
end