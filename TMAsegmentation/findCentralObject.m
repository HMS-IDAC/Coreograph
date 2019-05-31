function TMAmask = findCentralObject(bw)

TMAlabel = bwlabel(imresize(bw,0.25,'nearest'));
stats = regionprops(TMAlabel,'Centroid');

dist=[];
for iObject = 1: numel(stats)
    dist(iObject) = sqrt((stats(iObject).Centroid(2)-size(TMAlabel,1)/2)^2+(stats(iObject).Centroid(1)-size(TMAlabel,2)/2)^2);
end

% [minDistance, indexMin] = find(dist<0.6/2*sqrt(size(TMAlabel,1)*size(TMAlabel,2)));
% if isempty(indexMin)
    [minDistance, indexMin] = min(dist);    
% end
   
TMAmask = imfill(imclose(ismember(TMAlabel,indexMin),strel('disk',3)),'holes');