function [finalGrid] = gridFromCentroids(centroids, estCoreDiam, downsamplingFactor, varargin)

ip=inputParser;
ip.addRequired('centroids',@isnumeric);
ip.addOptional('showPlots',0,@isnumeric)
ip.parse(centroids,varargin{:});
showPlots=ip.Results.showPlots;

centroids = centroids.*(1/downsamplingFactor);

numCores=numel(centroids(:,1));
% estCoreDiam = 75;

%%plots the initial centroid positions on a black background
before = subplot(1,2,1);
before_pos = get(before,'Position');

maxX = round(max(centroids(:,1)))+200;
maxY = round(max(centroids(:,2)))+200;

if showPlots
    imagesub = zeros(maxY,maxX);
    imshow(imagesub,[])
    for iCore = 1:numCores
        hold on
        text(centroids(iCore,1),centroids(iCore,2),num2str(iCore),'Color','g')
    end
    before_pos(1) = 0.02;
    before_pos(3) = 0.47;
    set(before,'Position',before_pos);
end

%Picks out the core with coordinates closest to the middle of the TMA
%slide
midPoint = [maxX/2,maxY/2];
dist = sum((centroids(:,:)-midPoint) .^2,2);
[~, order] = sort(dist);
midCore = order(1);

%Retrieves all centroid distances from the core in the middle
dist = sum((centroids(:,:)-centroids(midCore,:)) .^2,2);
[sorted, order] = sort(dist);

%Goes through the nearest four neighboring cores and checks every pair of
%the neighboring cores for two of them that form a right angle to
%eachother
exit = 0;
for i = 2:4
    for j = 3:5
        n1 = (centroids(order(i),:)- centroids(midCore,:)) / norm(centroids(order(i),:) - centroids(midCore,:));
        n2 = (centroids(order(j),:) - centroids(midCore,:)) / norm(centroids(order(j),:) - centroids(midCore,:));
        if abs(90 - rad2deg(acos(dot(n1, n2)))) < 10
            closestCores = [order(i),order(j)];
            exit = 1;
            break
        end
    end
    if (exit)
        break
    end
end


%Save the centroids of the midCore and the two perpendicular neighbors
triangle=centroids([midCore, closestCores(1), closestCores(2)],:);

%Check which of the neigboring cores is above the midCore vertically and
%save it's index in the triangle as verti
if abs(triangle(1,2) - triangle(2,2)) > abs(triangle(1,2) - triangle(3,2))
    triangle(4,:) = [triangle(1,1) triangle(2,2)];
    verti = 2;
else
    triangle(4,:) = [triangle(1,1) triangle(3,2)];
    verti = 3;
end

%This portion is used in order to traverse down the column of the grid from
%the midCore in order to get a better estimate of what the vertical is for
%the TMA
distance = 0;
newPoint = triangle(verti,:); %Starts at the vertical closest core to midCore

%This ensures that the vector is always pointing in the downward direction
%by looking to see which has the bigger y value, the midcore or the
%vertical closest core found earlier.
if (triangle(verti,2) > triangle(1,2))
    n1 = (newPoint - triangle(1,:)) / norm(newPoint - triangle(1,:));
else
    n1 = (triangle(1,:) - newPoint) / norm(triangle(1,:) - newPoint);
end
prevPoint = newPoint;
order(1) = order(verti);

%This ensures that the vector that represents "down" is the best possible.
%It does this by traversing downwards looking for cores along the vector
%and adding their coordinates to the list and setting the final point to be
%an average of all of the points added to the list.
xVals = [];
yVals = [];
while (distance < estCoreDiam^2)
    
    %Translates point in "downward" direction
    translation = n1.*estCoreDiam;
    newPoint = prevPoint + translation;
    
    %Finds closest core and draws vector from previous point to closest
    %core
    dist = sum((centroids(:,:)-newPoint) .^2,2);
    [sorted, order] = sort(dist);
    distance = sorted(1);
    n2 = (centroids(order(1),:) - prevPoint) / norm(centroids(order(1),:) - prevPoint);
    
    %If the two vectors are within pi/8 radians of eachother, add the
    %closest cores coordinates to the list of coordinates
    if  (atan2(norm(det([n1;n2])),dot(n1,n2)) < pi/8)
        prevPoint = centroids(order(1),:);
        xVals = [xVals,centroids(order(1),1)];
        yVals = [yVals,centroids(order(1),2)];
        newPoint = [mean(xVals) mean(yVals)];
        if isequal(newPoint,centroids(midCore,:))
            continue
        end
        n1 = (newPoint - triangle(1,:)) / norm(newPoint - triangle(1,:));
        
        %Draws line on first plot to help visualize
        hold on
        line([newPoint(1),triangle(1,1)],[newPoint(2),triangle(1,2)],'Color', 'w');
    else
        break
    end
end

% Now I just use the previously found vector n2 as my "down/vertical"
% vector, previously [0 1]
vertical = (newPoint - triangle(1,:)) / norm(newPoint - triangle(1,:));

%Vector of remaining coordinates along with initial indices. Cores are
%removed once they are placed in the matrix, and removes midCore from the
%list
allCores = centroids.';
allCores(3,:) = 1:length(allCores);
remaining = allCores;
numPlaced = 1;

%Since there are sometimes disjoint sets on a single TMA, this first loop
%tries to fit the remaining cores to a grid until there are no more cores
%to place. This results in a list of grids that then need to be aligned
%together
grids={};
while numPlaced ~= 0 && length(remaining) > 0
    numRemaining = length(remaining);
    [remaining, placed] = fitToGrid(remaining, vertical, estCoreDiam);
    numPlaced = numRemaining - length(remaining);
    if (numPlaced > 0)
        grids{end+1} = placed;
    end
end



%This loop goes through all of the seperate grids and tries to align them
%and join them into one larger grid
while size(grids,2) > 1
    %Checks which one is the larger of the two grids and makes that the
    %main grid
    if (sum(size(grids{1})) > sum(size(grids{2})))
        main = grids{1};
        mainIndex=1;
        secondary = grids{2};
        secondaryIndex=2;
    else
        main = grids{2};
        mainIndex=2;
        secondary = grids{1};
        secondaryIndex=1;
    end
    
    numInMain = sum(main~=0,'all');
    numInSecondary = sum(secondary~=0,'all');
    mainCores = main(main~=0);
    
    %Calculates the average coordinate of the cores to find relative
    %positions
    mainCenter = sum(allCores(1:2,mainCores),2)/numInMain;
    secondaryCenter = sum(allCores(1:2,secondary(secondary~=0)),2)/numInSecondary;
    
    %Chooses which row and column of the secondary grid to use for
    %alignment base on where the secondary grid is relative to the main
    %grid
    row=size(secondary,1);
    column=size(secondary,2);
    toLeft = 1;
    above = 1;
    if (mainCenter(1) < secondaryCenter(1))
        column = 1;
        toLeft = 0;
    end
    if (mainCenter(2) < secondaryCenter(2))
        row = 1;
        above = 0;
    end
    
    secRow = size(secondary,1);
    secCol = size(secondary,2);
    while (secondary(row,column) == 0)
        if (row == 1 + (above)*(secRow-1))
            row = above*secRow+~above + ((2*~above)-1)*(toLeft*(secCol+1)+(2*~toLeft-1)*column);
            column = secCol*(toLeft) + ~toLeft;
        else
            row = row - (2*~above) + 1;
            column = column + (2*~toLeft) - 1;
        end
    end
    
    toAlign = allCores(:,secondary(row,column)); %Core used to align
    
    %Extends vertical line out from core to the edges of image to make the
    %vertical line
    translationFactor = (toAlign(2) - 1)/vertical(2);
    topPoint = toAlign(1:2).' - vertical*translationFactor;
    translationFactor = (maxY - toAlign(2))/vertical(2);
    bottomPoint = toAlign(1:2).' + vertical*translationFactor;
    
    %https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
    verticalDistance = sqrt(sum((bottomPoint-topPoint) .^ 2, 2));
    distances = zeros(1,numInMain);
    for i = 1:numInMain
        core = allCores(:,mainCores(i)).';
        area = abs((topPoint(2)-bottomPoint(2))*core(1) - ...
            (topPoint(1)-bottomPoint(1))*core(2) + ...
            topPoint(1)*bottomPoint(2) - bottomPoint(1)*topPoint(2));
        distances(i) = area/verticalDistance;
    end
    [sorted, order] = sort(distances);
    
    %Cores that are close to the line are considered further and used to
    %figure out which column of the main grid to place the core in
    candidates = mainCores(order(sorted < estCoreDiam*.35));
    if (~isempty(candidates))
        [~,cols] = find(ismember(main,candidates));
        newCol = mode(cols);
    elseif (toLeft)
        newCol = 0;
    else
        newCol = size(main,2)+1;
    end
    
    %Repeats process for the horizontal line
    horizontal = [-vertical(2), vertical(1)];
    translationFactor = (toAlign(1) - 1)/horizontal(1);
    leftPoint = toAlign(1:2).' - horizontal*translationFactor;
    translationFactor = (maxX - toAlign(1))/horizontal(1);
    rightPoint = toAlign(1:2).' + horizontal*translationFactor;
    horizontalDistance = sqrt(sum((rightPoint-leftPoint) .^ 2, 2));
    distances = zeros(1,numInMain);
    for i = 1:numInMain
        core = allCores(:,mainCores(i)).';
        area = abs((leftPoint(2)-rightPoint(2))*core(1) - ...
            (leftPoint(1)-rightPoint(1))*core(2) + ...
            leftPoint(1)*rightPoint(2) - rightPoint(1)*leftPoint(2));
        distances(i) = area/horizontalDistance;
    end
    [sorted, order] = sort(distances);
    candidates = mainCores(order(sorted < estCoreDiam*.35)); 
    if (~isempty(candidates))
        [rows,~] = find(ismember(main,candidates));
        newRow = mode(rows);
    elseif (above)
        newRow = 0;
    else
        newRow = size(main,1)+1;
    end
    
    %If the addition of the new grid is out of the boundarys of the main
    %grid it pads the matrix with zeros so that it can paste in the
    %secondary grid
    index = [newRow-(row-1),newCol-(column-1)];
    if (index(1) <= 0)
        main = [zeros(abs(index(1))+1,size(main,2)); main];
        index(1) = 1;
    elseif (index(1)+size(secondary,1)-1 > size(main,1))
        main = [main; zeros(index(1)+size(secondary,1)-1 - size(main,1),size(main,2))];
    end
    if (index(2) <= 0)
        main = [zeros(size(main,1),abs(index(2))+1) main];
        index(2) = 1;
    elseif (index(2)+size(secondary,2) > size(main,2))
        main = [main zeros(size(main,1),index(2)+size(secondary,2)-1 - size(main,2))];
    end
    
    
    %Pastes in the grid and removes secondary grid from list of grids
     for i = 1:numel(secondary)
         if (secondary(i) ~= 0)
             [row,col] = ind2sub(size(secondary),i);
             main(index(1)+row-1,index(2)+col-1) = secondary(i);
         end
     end
    
    
    %main(index(1):index(1)+(size(secondary,1)-1),index(2):index(2)+(size(secondary,2)-1)) = secondary;
    grids{mainIndex}=main;
    grids(:,secondaryIndex)=[]; 
end


if showPlots
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    after = subplot(1,2,2);
    imshow(imagesub,[])
    tmaGrid = grids{1};
    for row = 1:size(tmaGrid,1)
        for column = 1:size(tmaGrid,2)
            hold on
            coreNum = tmaGrid(row,column);
            if (coreNum == 0)
                continue
            end
            name = [Alphabet(row) num2str(column)];
            text(centroids(coreNum,1),centroids(coreNum,2),name,'Color','g')
        end
    end
    after_pos = get(after,'Position');
    set(after,'Position',[before_pos(1)+before_pos(3)+.02 before_pos(2) .47 before_pos(4)]);
end

finalGrid = grids{1};
end


function [remaining, matches] = fitToGrid(corePositions, vertical, estCoreDiam)
maxX = round(max(corePositions(1,:)))+200;
maxY = round(max(corePositions(2,:)))+200;

midPoint = [maxX/2;maxY/2];
dist = sum((corePositions(1:2,:)-midPoint) .^2,1);
[~, order] = sort(dist);
midCore = order(1);

%Initializes the data structures
index = [1 1];
stack = [];
remaining = corePositions;
prev = [remaining(:,midCore).' index 0];
grid = [];
grid(1,1,:) = remaining(1:2,midCore);
matches = [remaining(3,midCore)]; %Puts midcore at (1,1) to start
remaining(:,midCore) = [];
prevIndex = prev(4:5);
dist = sqrt(sum((remaining(1:2,:)-prev(1:2).') .^ 2, 1));
[sorted, order] = sort(dist);
remain = length(order);
stack = [stack;[remaining(:,order).' repmat(prevIndex,remain,1) sorted.']];
stack = sortrows(stack, 6);


reset = 0;

%Goes until there are no more cores to place on the grid
while (~isempty(stack))
    %Pops the shortest found distance of the stack
    next = stack(1,:);
    stack(1,:) = [];
    
    %If the distance is too great stop trying to place them
    if (next(6) > 1.5*estCoreDiam)
        break
    end
    
    %Grabs both the previous index and creates a vector of the previous
    %core
    prevIndex = next(4:5);
    prev = [grid(prevIndex(1),prevIndex(2),1) grid(prevIndex(1),prevIndex(2),2) matches(prevIndex(1),prevIndex(2)) prevIndex];
    
    %Calculates the normalized vector from previous to next core
    n1 = (next(1:2) - prev(1:2)) / norm(next(1:2) - prev(1:2));
    if (next(3) == 80)
        %disp("helo") Used for debugging certain cores
    end
    
    %Checks whether it is within range of a 90 degree position
    if (abs(90 - rad2deg(acos(dot(n1, vertical)))) < 25 || 270 - rad2deg(acos(dot(n1, vertical))) < 25)
        if (next(1) - prev(1)) > 0
            index = [prevIndex(1) prevIndex(2)+1];
        else
            index = [prevIndex(1) prevIndex(2)-1];
        end
    elseif (rad2deg(acos(dot(n1, vertical))) < 25 || 180 - rad2deg(acos(dot(n1, vertical))) < 25)
        if (next(2) - prev(2)) > 0
            index = [prevIndex(1)+1 prevIndex(2)];
        else
            index = [prevIndex(1)-1 prevIndex(2)];
        end
        %If it is not close to one of the 90 degree positions draw a line
        %opposite of the vector from prev to next and see if a core has
        %already been placed in the opposite position (allows diagonal
        %placement).
    else
        %Finds core most likely in opposite position
        opposite = prev(1:2) - (n1*estCoreDiam);
        dist = sqrt(sum((corePositions(1:2,:)-opposite(1:2).') .^ 2, 1));
        [sorted, order] = sort(dist);
        if (sorted(1) < estCoreDiam*.5)
            oppositeCore = order(1);
        else
            continue
        end
        
        %If the opposite core has already been placed make vector from
        %previous core to the opposite core
        oppositeIndex = find(matches == oppositeCore);
        if ~isempty(oppositeIndex)
            [oppositeIndex(1), oppositeIndex(2)] = ind2sub(size(matches),oppositeIndex);
            oppositeCoreCoords = [corePositions(1:2,oppositeCore)].';
            prevVec = (oppositeCoreCoords(1:2) - prev(1:2)) / norm(oppositeCoreCoords(1:2) - prev(1:2));
        else
            continue
        end
        
        %If the angle between the preVec and n1 are close to 180 then place
        %the next core in the position opposite of the opposite core
        if (180 - rad2deg(acos(dot(n1, prevVec))) < 25)
            changeInIndex = prevIndex - oppositeIndex;
            index = prevIndex + changeInIndex;
        else
            continue
        end
    end
    
    %If the new index is in a new row/column pad the matches array with
    %zeros in that new row or column and adjust all the previous indices
    if (index(1) == 0)
        matches = [zeros(1,size(matches,2)); matches];
        grid = [zeros(1,size(grid,2),2);grid];
        stack(:,4) = stack(:,4) + 1;
        index(1) = 1;
    elseif (index(2)== 0)
        matches = [zeros(size(matches,1),1) matches];
        grid = [zeros(size(grid,1),1,2) grid];
        stack(:,5) = stack(:,5) + 1;
        index(2) = 1;
    elseif (index(1) > size(matches,1))
        matches = [matches; zeros(1,size(matches,2))];
        grid = [grid;zeros(1,size(grid,2),2)];
    elseif (index(2) > size(matches,2))
        matches = [matches zeros(size(matches,1),1)];
        grid = [grid zeros(size(grid,1),1,2)];
    elseif (matches(index(1),index(2)) ~= 0)
        continue
    end
    
    %Place the core in the grid and remove it from the stack and list of
    %remaining cores
    next(4:5) = index;
    matches(index(1),index(2)) = next(3);
    grid(index(1),index(2),:) = next(1:2);
    prev = next;
    remaining(:,remaining(3,:) == next(3)) = [];
    stack(stack(:,3) == next(3),:) = [];
    
    
    %Calculate the distances to the remaining cores and place them on the
    %stack for next iteration
    prevIndex = prev(4:5);
    dist = sqrt(sum((remaining(1:2,:)-prev(1:2).') .^ 2, 1));
    [sorted, order] = sort(dist);
    remain = length(order);
    stack = [stack;[remaining(:,order(1:remain)).' repmat(prevIndex,remain,1) sorted(1:remain).']];
    stack = sortrows(stack, 6);
end
end
