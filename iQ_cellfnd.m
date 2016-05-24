function [DIC_bd,mask_DIC_grid,DIC_coor]=iQ_cellfnd(DIC,dil_sz,savefig)
%% THIS FUNCTION FINDOUT THE POSITIONS AND BOUNDARIES OF EACH CELL AND GENERATE THE CORRECT MASK FOR FLUORESCENCE IMAGES
% INPUTS:
% DIC: DIC image of target cell
% dil_sz: dilate size of the cell boundary
% savefig: input 1 if want to save figure into ppt, otherwise left it or input 0.

% OUTPUTS:
% DIC_bd: an overlay image of DIC cell with calculated DIC cell boundary  
% mask_DIC_grid: two column x,y coordinates inside each cell mask.
% DIC_coor: two column x, y coordinates of each cell mask boundary.

%% save ppt or not
if nargin<3;    savefig=0; end;
if savefig==1; pptfilname=('DIC.ppt'); end;
currentdic=pwd;

%% GENERATE DIC MASK FOR EACH CELL
I1=DIC;        % load the image
background = imopen(I1,strel('disk',15));  % calculate the bkg using strel function

I2 = I1- background; % remove the background so that you can have a even intensity image.
int_th=mean(mean(I2))+3*mean(std(double(I2)));          
I2(I2<int_th)=0;     % remove too weak intinsity stuff                          
                      
% finding the edge using sobel alogrithm which find the maximum gradients                               
I3 = edge(I2,'sobel');
se1 = strel('disk',1);     
I3=imdilate(I3,se1);       % dilate the edge by 1 pixel to make sure the cell has a closed  boundary to fill
I3 = imfill(I3, 'holes');
I3=imerode(I3,se1);        % once the cell was filled, erode one pixel to get the correct size of the cell
close all

I4=imclearborder(I3, 4); % remove all objects connected to the border
[~,L] = bwboundaries(I4,'noholes');
stats = regionprops(L,'Area');
statsarea=cell2mat(struct2cell(stats));
hist(statsarea,50);      % generate a area histogram
xlabel('Area'); ylabel('Counts');

n=1;                     % just iterations to find out the best minimum cell area.
    while 1
    th_sz=input('min and max cell area [min, max]: ');
    th_szmin=th_sz(1); th_szmax=th_sz(2);
    I5min=bwareaopen(I4,th_szmin);
    I5max=bwareaopen(I4,th_szmax);
    I5=I5min-I5max;
        if n == 1;
            figure; subplot(121); imshow(I1,[]); title('original image');
        end
            subplot(122); imshow(I5,[]); title('too small object removed image');
    b=input('Want to Change Size Threshold ?(Y/N ) ','s');
    B=upper(b);
        if B=='Y'; n=n+1; end;
        if B=='N', break; end;
    end

% Dilate each mask so that the mask can cover the whole cell accordingly.
cc = bwconncomp(I5, 4); %Identify Objects in the Image
I6=false(size(I5));

CellDilMask=cell(cc.NumObjects,1);
CellDilPrim=cell(cc.NumObjects,1);
for i=1:cc.NumObjects
    cells= false(size(I5));
    cells(cc.PixelIdxList{i}) = true;
    se = strel('disk',dil_sz);
    cells_dil = imdilate(cells,se);
        I6=I6+cells_dil;
    cells_dil_perim=bwboundaries(cells_dil,'noholes');
    cells_dil_perim=cell2mat(cells_dil_perim);
    CellDilPrim{i}=cells_dil_perim;
    cells_dil_mask=find(cells_dil~=0);
    CellDilMask{i}=cells_dil_mask;
end

% choose cells of interest by the user
close all
figure; imshow(DIC,[]);
hold on
for b=1:cc.NumObjects
    plot(CellDilPrim{b}(:,2),CellDilPrim{b}(:,1));
    cellidx=cc.PixelIdxList{b}(1);
    [cell_y,cell_x]=ind2sub(size(I6),cellidx);
    metric_string = int2str(b);
    text(cell_x-15,cell_y-5,metric_string,'Color','y',...
       'FontSize',12,'FontWeight','bold');
end
hold off;
    c=input('Want to pick up mask yourself ?(Y/N ) ','s');
    C=upper(c);
    if C=='Y'; cellofinterest=input('Input cell number index= ');
        nocellint=length(cellofinterest);
        CellDilPrim_choose=CellDilPrim(cellofinterest);
        I7=false(size(I6));
        figure
        imshow(I7,[]); hold on; 
        for p=1:nocellint;
            plot(CellDilPrim_choose{p}(:,2),CellDilPrim_choose{p}(:,1));
            cell_y = CellDilPrim_choose{p}(1,1);
            cell_x = CellDilPrim_choose{p}(1,2);
            metric_string = int2str(p);
            text(cell_x-15,cell_y-5,metric_string,'Color','y',...
            'FontSize',12,'FontWeight','bold');
        end
        hold off;
    end;
    if C=='N', nocellint= cc.NumObjects;
       CellDilPrim_choose=CellDilPrim; 
    end;

% generate DIC image with boundary
    DIC_coor=fliplr(cell2mat(CellDilPrim_choose));
    DIC_bd=DIC;
    DIC_bd(sub2ind(size(DIC_bd),DIC_coor(:,2),DIC_coor(:,1)))=max(max(DIC_bd))+0.3*(max(max(DIC_bd))-min(min(DIC_bd)));
    [row col]=find(I6>0); 
    mask_DIC_grid=[col row];

%%  Plot out important images 
    close all    
    % plot all the mask in one figure 
    fh1=figure;                      
    subplot(231); imshow(I1,[]); title('original image');
    subplot(232); imshow(I2,[]); title('bkg removed image');
    subplot(233); imshow(I3,[]); title('threshold image');
    subplot(234); imshow(I4,[]); title('cleared border image');
    subplot(235); imshow(I5,[]); title('too small object removed image')
    uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
    
    fh2=figure; 
        imshow(DIC,[]); hold on; title('DIC image with boundary')
        for p=1:nocellint;
            plot(CellDilPrim_choose{p}(:,2),CellDilPrim_choose{p}(:,1));
            cell_y = CellDilPrim_choose{p}(1,1);
            cell_x = CellDilPrim_choose{p}(1,2);
            metric_string = int2str(p);
            text(cell_x-15,cell_y-5,metric_string,'Color','y',...
            'FontSize',12,'FontWeight','bold');
        end
        hold off;
    uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
    
    fh3=figure;
       ha = tight_subplot(ceil(nocellint/5),5,[.01 .01],[.1 .01],[.01 .01]);
       for q = 1:nocellint;
           axes(ha(q));
           imshow(I1,[]);
           cen_y = mean(CellDilPrim_choose{q}(:,1));
           cen_x = mean(CellDilPrim_choose{q}(:,2));
           axis([cen_x-25,cen_x+25,cen_y-25,cen_y+25]);
           metric_string = int2str(q);
           text(cen_x-15,cen_y-15,metric_string,'Color','y',...
            'FontSize',8,'FontWeight','bold');
       end
           aaa = (nocellint/5); bbb = ceil(nocellint/5);
           if aaa ~= bbb;
              exsubran=(nocellint+1:1:ceil(nocellint/5)*5);
              delete(subplot(ceil(nocellint/5),5,exsubran));
           end
      uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
      
      if savefig==1;
      saveppt2(pptfilname,'figure',fh1,'scale',true,'stretch',false);    
      saveppt2(pptfilname,'figure',[fh2,fh3],'scale',true,'stretch',false);
      end
end