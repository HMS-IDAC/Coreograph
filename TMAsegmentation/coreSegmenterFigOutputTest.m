function TMAmask = coreSegmenterFigOutputTest(DAPI,varargin)

ip = inputParser;
ip.addParamValue('activeContours','true',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('split','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('initialmask',[],@(x)(numel(x)>0));
ip.addParamValue('preBlur',[],@(x)(isnumeric(x) & numel(x)>0));
ip.parse(varargin{:});          
p = ip.Results;  

% nucGF=stdfilt(imresize(DAPI,0.5),ones(3,3));
nucGF=imgaussfilt(imresize(DAPI,0.5),3);
buffer = round(0.2*size(nucGF,1));

nuclearMask = imgaussfilt(nucGF,p.preBlur)>thresholdOtsu(nucGF);
stats=regionprops(nuclearMask,'Area','Solidity','Eccentricity');
[M,I]=max(cat(1,stats.Area));

if (stats(I).Solidity<0.9 || stats(I).Area/size(nuclearMask,1)/size(nuclearMask,2) <0.2) || isequal(p.activeContours,'true')
    if isempty(p.initialmask)
        p.initialmask = zeros(size(nucGF));
        p.initialmask(buffer:end-buffer,buffer:end-buffer) = 1;
    end
    initialmask = imresize(p.initialmask,size(nucGF));
    nuclearMask = activecontour(nucGF,initialmask,300,'Chan-Vese','SmoothFactor',2);
    stats=regionprops(nuclearMask,'Area','Solidity','Eccentricity');
    [M,I]=max(cat(1,stats.Area));
end
    
%% watershed segmentation
if isequal(p.split,'true') && ((sum(sum(nuclearMask))/size(nucGF,1)/size(nucGF,2)>0.4) || stats(I).Eccentricity>0.8) %assume the object is the largest one
%     nMaskDist =imgaussfilt3(-bwdist(~nuclearMask),2*p.preBlur);
    nMaskDist =imhmin(-bwdist(~nuclearMask),2*p.preBlur);
    cytograd= imimposemin(nMaskDist,imerode(~nuclearMask,strel('disk',3))| imregionalmin(nMaskDist));
    TMAmask=watershed(cytograd);
    TMAmask = nuclearMask.*(TMAmask>0);
    
else
    TMAmask = nuclearMask;
end
%  
TMAmask = bwareaopen(TMAmask,round(size(TMAmask,1)*size(TMAmask,2)*0.005));
TMAlabel = bwlabel(TMAmask);
stats= regionprops(TMAlabel);

%% choose the object closest to center of image
dist=[];
for iObject = 1: numel(stats)
    dist(iObject) = sqrt((stats(iObject).Centroid(1)-size(nucGF,1)/2)^2+(stats(iObject).Centroid(2)-size(nucGF,2)/2)^2);
end

[minDistance, indexMin] = find(dist<0.6/2*sqrt(size(TMAlabel,1)*size(TMAlabel,2)));
if isempty(indexMin)
    [minDistance, indexMin] = min(dist);  
else 
    [minDistance, indexMin] = sort(dist);  
end
   
TMAmask = imfill(imclose(ismember(TMAlabel,indexMin(1)),strel('disk',3)),'holes');
%  imshowpair(bwperim(imresize(TMAmask,size(DAPI))),imadjust(DAPI))


