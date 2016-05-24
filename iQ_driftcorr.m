function [ DFT_pos,DFT_dr,DFT_dr_ave ] = iQ_driftcorr( mv,threshold,mask_sz,IM_info,refcoor,Z_cal_curve )
%% THIS FUNCTION ESTIMATE THE AVERAGE DRIFT RELATIVE TO THE FIRST FRAME OF THE MOVIE FOR ASSIGNED SPOTS
% INPUT:
% mv: movie folder
% threshold: minimum number of sigma needed for a peak picking, suggested value is 4. 
% mask_sz: size of mask that cover each spots, recommend size:13
% IM_info: contains pixel size information for x and y,ex. IM_info=[135.8 176.6]
% refcoor: you can manually input the coordinates of marker of interest in [x1,y1;x2,y2;x3,y3] format.
% Z_cal_curve: Z calibration curve generated from iQ_ZCal.m

% note: if you are only dealing w/ 2D data, you can just input the first
% three inputs, the program will recognize that you are doing 2D data analysis. 

% OUTPUTS:
% DFT_pos: center positions of picked spots through the movie, each row
% contains following information (16 colums in total): 
% [mu_x(px) mu_y(px) mu_x(nm) mu_y(nm) mu_z(nm) 
% sigma_x(px) sigma_y(px) sigma_x(nm) sigma_y(nm) ellip PSFW_H(px) 
% n amp background Int std_b] by 2d gaussian fit;
% DFT_dr: drift estimation of each picked spot
% DFT_dr_ave: average drift estimation of picked spots

%% Top of the routine
tic
cd(mv); D=dir;
num_frame=length(D)-2;
im1=double(imread(D(3).name));
if isempty(refcoor)
%%   FIND OUT THE POSSIBLE SPOTS FOR DRIFT CORRECTION
%  check out 6 individual images of the movie, the spacing between each
%  image is equal to (frame number of movie/6). The purpose of this
%  process is to determine which spots are ideal for drift correction.
plotlist=(3:floor(num_frame/6):num_frame-2)';
    figure
    for m=1:length(plotlist)
    plotnum=num2str(plotlist(m,1)-2);
    plotname=strcat('frame no. ',plotnum);
    im=imread(D(plotlist(m,1)).name);
    subplot(2,3,m);
    imshow(im,[]);
    title(plotname);
    end
    
% After going through 6 frames, ask for procedding to pick up spots.    
    PS=input('pick up the spots (Y)? ','s');
    PS=upper(PS);
    if isempty(PS), PS='Y';end;
        if PS=='Y', disp('Continue')
    end
        
%% PICK THE SPOTS FOR DRIFT CORRECTION
% this process decides the spots for drift correction by clicking on the
% spots on the image averaged over first 20 frames.

% load and average the first 20 frames to determine the calibration spots
im_pick=double(imread(D(3).name));
for k=2:20; 
 image=double(imread(D(k+2).name));                     
 im_pick=im_pick+image;
end
im_pick=im_pick/20;
% show the image 
figure;  imshow(im_pick,[])                              
% pick up spots for drift correction by clicking on the spots
spot_coor=putp;
pick_hand=spot_coor(:,1); pick_hand(:,2)=spot_coor(:,2);
else
    pick_hand = refcoor;
end

% find out the local maximum spots of frame No.1
 [~,peak_pick,~,~]=iQ_pkfnd(im1,threshold,mask_sz);             
 
% find out the exact coordinates of assigned spots using function iQ_nearestPoints
 pick= iQ_nearestPoints(pick_hand,peak_pick);               
 pick=pick';  % transpose pick to make it has [x1 y1 x2 y2 x3 y3] format intead of [x1 y1;x2 y2;x3 y3] to meet the input format of iQ_gaussianfit
 pick=cell2mat(pick);
 
DFT_CORR=cell(num_frame,1); 

%% FIND OUT PEAK POSITIONS OF SPOTS EVERY FRAME THROUGH THE MV
% load every frame for the calibration curve 
for pki=1:num_frame;
    display(['Frame no. ',int2str(pki)])
    im=double(imread(D(pki+2).name));
    
% fit the background and remove the background from original image    
     bk_fit = imopen(im,strel('disk',15));               
     im_bk=im-bk_fit;                                    

% get xysxsy == [mu_x(px) mu_y(px) mu_x(nm) mu_y(nm) mu_z(nm) 
%                sigma_x(px) sigma_y(px) sigma_x(nm) sigma_y(nm) ellip PSFW_H(px)
%                n amp background Int std_b] by 2d gaussian fit;
     if nargin == 6;
     xysxsy = iQ_gaussianfit( pick,im_bk,mask_sz,pki,IM_info,Z_cal_curve );
     else
     xysxsy = iQ_gaussianfit( pick,im_bk,mask_sz,pki,IM_info);
     end
     xysxsy=cell2mat(xysxsy);
     DFT_CORR{pki,1}=xysxsy;     
      
end

%% GET THE DRIFT DISTANCE COMPARING TO THE FIRST FRAME
 DFT_pos=cell2mat(DFT_CORR);
 DFT_dr=DFT_pos;
 
 [~,c]=size(DFT_pos);
 % calculate the drift by substrating the first row of each spot. 
 for k=1:c/16;
     DFT_dr(:,16*(k-1)+1)=DFT_pos(:,16*(k-1)+1)-DFT_pos(1,16*(k-1)+1);
     DFT_dr(:,16*(k-1)+2)=DFT_pos(:,16*(k-1)+2)-DFT_pos(1,16*(k-1)+2);
     DFT_dr(:,16*(k-1)+3)=DFT_pos(:,16*(k-1)+3)-DFT_pos(1,16*(k-1)+3); 
     DFT_dr(:,16*(k-1)+4)=DFT_pos(:,16*(k-1)+4)-DFT_pos(1,16*(k-1)+4);
     DFT_dr(:,16*(k-1)+5)=DFT_pos(:,16*(k-1)+5)-DFT_pos(1,16*(k-1)+5); 
     
 end
 
%% AVERAGE OVER SELECTED SPOTS TO GET  AVERAGE DRIFT CORRECTION
% only one set of reference point case 
 if c==16;
 DFT_dr_ave=DFT_dr(:,1:5);
 end
 
% multiple sets of reference points case 
     DFT_dr_ave=DFT_dr(:,1:5);
 for o=2:c/16;
     DFT=DFT_dr(:,16*(o-1)+1:16*(o-1)+5);
     DFT_dr_ave=DFT_dr_ave+DFT;
 end
 DFT_dr_ave=DFT_dr_ave/(c/16);
 toc
%% nested function iQ_pkfnd
%% nested function iQ_nearestPoints
%% nested function iQ_gaussianfit
%% nested function putp
end

