function [NPs]=iQ_nearestPoints(peak,target)
% Purpose: try to find out the nearest points related to peaks in target
% point sets and output the points as NPs
% INPUT: peak, target
%        peak: the points you pick by your hand
%        targets: target peaks that contain many points
% OUTPUT:NP
%        NPs: nearestPoints, one column cell data.
x=target(:,1);
y=target(:,2);
[row ~]=size(peak); NPs=cell(row,1);
for npi =1:row;
d = sqrt((target(:,1) - peak(npi,1) ).^2 + (target(:,2) - peak(npi,2) ).^2 );

[~, ixSort] = sort(d); %Sort distances

nearestpoint=[x(ixSort(1)) y(ixSort(1))] ;
NPs{npi,1}=nearestpoint;
end
end
