clc;clear all;close;

moviefile(1) = cellstr('D:\Research\2012 fall\ZXC_Movie20120824\1_p10_rod_800uM_AR_60mM_H2O2_Selector3.tif');
moviefile(2) = cellstr('D:\Research\2012 fall\ZXC_Movie20120824\1_p10_rod_800uM_AR_60mM_H2O2_Selector4.tif');
% moviefile(3) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector3.tif');
% moviefile(4) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector4.tif');
% moviefile(5) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector5.tif');
% moviefile(6) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector6.tif');
% moviefile(7) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector7.tif');
% moviefile(8) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector8.tif');
% moviefile(9) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector9.tif');
% moviefile(10) = cellstr('D:\Research\2012 fall\11012012gnr150nmRZ\movies\20121101 150nm rz_Selector10.tif');
selectorlength=10000;
number_avg=5;
height=200;
width=400;
a=zeros(height,width);
avg=zeros(height,width);
for i=1:length(moviefile)
moviename = char( moviefile(i) );
for j=1:number_avg
 a(:,:) = imread(moviename,j);
 a=double(a);
 avg=avg+a;
end
end