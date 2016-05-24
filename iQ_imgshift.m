function [Red_DIC, mask_Red_grid,Rshf,DCA,RCA ] = iQ_imgshift( DIC,Red,DIC_coor,mask_sz,IM_info,refcoor,Z_cal_curve)
%% THIS FUNCTION IS TO SHIFT THE MASK COORDINATES TO GENERATE PROPER MASK FOR FLUORESCENCE IMAGES
%   INPUTS:
%       DIC: DIC image
%       Red: Red image
%       DIC_coor: two column data for the boundary coordinates of a cell
%       mask_sz: mask size for PSF fitting, typically use 13
%       IM_info: pixels information, two element array: [pixel x, pixel y];
%       refcoor: two column data specify reference point coor, if you want to skip clicking on the
%       image to specify the coor, input the coordinates here in [x1,y1;x2,y2;x3,y3] format.
%       Z_cal_curve: z calibration curve
%   OUTPUTS:
%       Red_DIC: overlay image of Red image with cell mask boundary
%       mask_Red_grid: two column data for all coordinates inside each cell mask
%       Rshf: unit: pixel
%       Rshf=[Red_xshf,Red_yshf,mean(DIC_Au_x-Red_Au_x)+dx_R,mean(DIC_Au_y-Red_Au_y)+dy_R];
%       4 element array, Rshf(1) and (2) give the rounded image shift x and y used for cell mask shifting
%                        Rshf(3) and (4) give the exact image shift x,y 
%       DCA: unit: pixel
%       DCA(k,1)=xbar; DCA(k,2)=ybar;
%       2 column matrix, column 1 and 2 give the centroid position of kth reference point 
%                        
%       RCA: unit: pixel
%       RCA=[mu_x,mu_y,sigma_x,sigma_y,PSFW_H, mu_z];
%       5 column matrix, PSF fitting result of reference point, column 1 to
%       5 give ceter x,y, sigmax,sigmay,PSF width-height and correspoding z
%       position. Note: here mu_y and sigma_y have been corrected for
%       cylindrical distortion.
%% TOP of the routine
close all;
DIC = double(DIC); DIC = DIC-min(min(DIC));
bkDIC = imopen(DIC,strel('disk',15));
DIC=DIC-bkDIC;

Red = double(Red); Red = Red-min(min(Red));
bkRed = imopen(Red,strel('disk',15));
Red=Red-bkRed;

if isempty(refcoor) ;

h1= figure;   imshow(DIC,[]); imcontrast;
h2= figure;   imshow(Red,[]); imcontrast; title('Red image');

% choose reference spots
a=1; b=1;
    while 1
        q1=input('Have enough pick spots ?(Y/N ) ','s');
        Q1=upper(q1);
        if Q1=='Y'; hold off; break; end;
        if Q1== 'N'; a = a+1; end;
        
        figure(h1); hold on;
        ax=axis; xspacing=(ax(1,2)-ax(1,1))/100; yspacing=(ax(1,4)-ax(1,3))/100;
        [pkx,pky]=ginput(1);
        plot(pkx,pky,'rx'); number=num2str(b);text(pkx+xspacing,pky+yspacing,number,'color','w')
        
        figure(h2);  hold on;
        plot(pkx,pky,'yo'); number=num2str(b);text(pkx+xspacing,pky+yspacing,number,'color','w')
        
        pkall(b,:)=[pkx,pky];
        b=b+1;
    end
pkallrnd=round(pkall); pkallrnd=sortrows(pkallrnd,1);
else
    pkallrnd=sortrows(refcoor);
end

% only look picked spots and do the sequential analysis
%DIC case
nDIC = false(size(DIC));
nDIC(sub2ind(size(nDIC),pkallrnd(:,2),pkallrnd(:,1)))=1;
nDIC1=imfilter(nDIC,ones(5)); 
nDIC=nDIC1.*DIC; imtool(nDIC,[]); 

    L_au = logical(nDIC1);
    s_au = regionprops(L_au,'PixelIdxList', 'PixelList');
    DCA=zeros(numel(s_au),2); 
    for k=1:numel(s_au);
        idx= s_au(k).PixelIdxList;
        pixel_values = double(nDIC(idx));
        sum_pixel_values = sum(pixel_values);
        x = s_au(k).PixelList(:, 1);   y = s_au(k).PixelList(:, 2);
        xbar = sum(x .* pixel_values) / sum_pixel_values;
        ybar = sum(y .* pixel_values) / sum_pixel_values;       
        DCA(k,1)=xbar; DCA(k,2)=ybar;
    end
    
DIC_Au_x=DCA(:,1); DIC_Au_y=DCA(:,2);
    
% find out the centroid,area, major/minor axis information for Red
%Red case
nRed = false(size(Red));
nRed(sub2ind(size(nRed),pkallrnd(:,2),pkallrnd(:,1)))=1;
nRed=imfilter(nRed,ones(21));
nRed=nRed.*Red;

pixel_x=IM_info(1); pixel_y=IM_info(2);  corr_cylin=pixel_x/pixel_y;
[row_idxau,~]=size(pkallrnd);
sz_red=ceil(mask_sz/2);
RCA=cell(row_idxau,1);
    for jj=1:row_idxau; 
        [xi,yi]=meshgrid(pkallrnd(jj,1)-sz_red:pkallrnd(jj,1)+sz_red,pkallrnd(jj,2)-sz_red:pkallrnd(jj,2)+sz_red);
        zi=nRed(pkallrnd(jj,2)-sz_red:pkallrnd(jj,2)+sz_red,pkallrnd(jj,1)-sz_red:pkallrnd(jj,1)+sz_red);
        results = iQ_autoGaussianSurfML(xi,yi,zi);
        mu_x=results.x0; mu_y=results.y0/corr_cylin;
        sigma_x=results.sigmax; sigma_y=results.sigmay/corr_cylin;
        PSFW_H=sigma_x-sigma_y;
        mu_z=0;
        if nargin==7;
        x=Z_cal_curve(:,1);                                     
        y=Z_cal_curve(:,2);
        % calulate z from the width-height of psf
        mu_z=interp1(y,x,PSFW_H);
        end
        
        RCA{jj}=[mu_x,mu_y,sigma_x,sigma_y,PSFW_H, mu_z];
    end
RCA=cell2mat(RCA);
Red_Au_x=RCA(:,1); 
Red_Au_y=RCA(:,2)*corr_cylin;  % here multiply corr_cylin is due to we want to keep y pixel size relative to x the same in DIC case

% correct image coordinates with Au NP coordinates
% calculate the average deltax and deltay for correction
dx_R = 0; dy_R = 0;
    while 1
    Red_xshf=round(mean(DIC_Au_x-Red_Au_x))+dx_R; Red_yshf=round(mean(DIC_Au_y-Red_Au_y))+dy_R;
    
    Red_coor_x=DIC_coor(:,1)-Red_xshf;    % correct for x
    Red_coor_y=DIC_coor(:,2)-Red_yshf;    % correct for y
    
% if new boundary coordinates after correction is outside the image range, trim it to the boundary of the image
    [max_y max_x]=size(DIC);
    Red_coor_x(Red_coor_x<1)=1; Red_coor_x(Red_coor_x>max_x)=max_x;
    Red_coor_y(Red_coor_y<1)=1; Red_coor_y(Red_coor_y>max_y)=max_y;
    Red_coor=[Red_coor_x,Red_coor_y];
    
% use new boundary coordinates to create new mask grid
    mask_Red=zeros(size(DIC));
    mask_Red(sub2ind(size(mask_Red),Red_coor(:,2),Red_coor(:,1)))=1;
    mask_Red=imfill(mask_Red,'holes');
    mask_Red=imclearborder(mask_Red, 4); % remove the mask attach the image boundary

% get the coordinates of every point of the mask grid
    [row col]=find(mask_Red>0); 
    mask_Red_grid=[col row];
% get the boundary coordinates of new mask that boundary masks have been removed
    mask_Red_bd=bwboundaries(mask_Red,'noholes');
    mask_Red_bd=cell2mat(mask_Red_bd);
    
% generate the overlay image of Red image with corrected DIC boundary
    Red_DIC=Red;
    Red_DIC(sub2ind(size(Red_DIC),mask_Red_bd(:,1),mask_Red_bd(:,2)))=max(max(Red))*1.5;
    imshow(Red_DIC,[]); imcontrast;
    b=input('need to fine tune the position ?(Y/N ) ','s');
    B=upper(b); 
    if B=='Y'; 
       xyfine_R=input('[\Deltax,\Deltay] (must be integer)=');
       dx_R = xyfine_R(1); dy_R = xyfine_R(2);
        n=n+1; 
    end;    
    if B=='N', break; end;
    end
    
    Rshf=[Red_xshf,Red_yshf,mean(DIC_Au_x-Red_Au_x)+dx_R,mean(DIC_Au_y-Red_Au_y)+dy_R];

end


