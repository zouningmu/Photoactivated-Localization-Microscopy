function [RI]=iQ_ReconImage(Axis_coor,rangexy,mask_size,sigma,precision)
% this function is for reconstruct PALM image with specified Gaussian function.
% INPUTS:
% Axis_coor: files contains x,y coordinates
% rangexy: two element array, specify the range of reconstruted image
% mask_size: deifine the mask size of the gaussian filter.
% sigma: define the sigma of the gaussian. 
% precision: how many devisions you want to have in the original pixel
% OUTPUS:
% RI: recontructed image

%% TOP of the routine
tic
if nargin < 5, precision = 10; end             %  precision: this input specifies the precision you want to have, if
                                               %  you want to have precision down to one decimal point, use 10. for
                                               %  the PALM image purpose 10 is ideal.
if nargin < 4, sigma = 1; end
if nargin < 3, mask_size = 11; end


rangx=rangexy(1); rangy=rangexy(2);

x=Axis_coor(:,1);  max_x=rangx;
y=Axis_coor(:,2);  max_y=rangy;

minx=min(x); miny=min(y);
if minx<0 | miny<0;
    display('Error: Axis_coor contains negative values !! Recheck your Axis_coor file')
    xy=[x,y];
    xy(xy<0)=0;
    xy(any(xy==0,2),:)=[];
    x=xy(:,1); y=xy(:,2);
end
% mi = min(min(x),min(y)); 
% x = x - mi;             % subtract minimum to get rid of negative numbers
% y = y - mi;
x = round(x*precision);    
y = round(y*precision);    

%% CONSTRUCT BLANK HIGH RESOLUTION IMAGE
img = zeros(round(max_y*precision),round(max_x*precision));    %# image will be square, doesn't have to be
x = uint32(x);
y = uint32(y);

%% MARK THE POSITION WITH TARGET XY DATA
ind = sub2ind(size(img),y,x);    %# use x,y or reverse arrays to flip image
img(ind) = 1;    %# could set intensity or RGB values in a loop instead

%% CONVOLUTE THE TARGET XY WITH PSF FUNCTION
RI=imfilter(img,fspecial('gaussian',mask_size,sigma));
%imshow(RI,[]);
toc
end
