%This file contains the function that takes in a single field from a lif
%file, finds the cilia in the field based on user-defined
%parameters/thresholds, and quantifies the signal of the
%protein-of-interest, such as Gli2, which is normally located in a specific
%part of the cilium
function [true_cilia_indices,possb_cilia_info,poi_intensities,poi_areas]=quant_Gli2_for_field(field,PARAMS)

cilia_channel=field{PARAMS.ciliaNum,1}; %get max-z projection values for the cilia staining

poi_channel=field{PARAMS.targetNum,1}; %get max-z projection values for the protein-of-interest (POI) staining

%get rid of salt and pepper noise in the first channel and then create a
%black and white mask by thresholding
M=(medfilt2(cilia_channel,[3 3])>PARAMS.Threshold); %identify pixels with real signal and store them in M

%save the cilia mask as tif file for troubleshooting
imwrite(M, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_mask.tif']);

%identify and number each object that could be a cilium in the mask.
%Objects are identified based on adjacent pixels that show some signal in
%the mask.
L=bwlabel(M);

%define the features/thresholds of objects that are cilia-like
MinArea=PARAMS.MinArea;
MaxArea=PARAMS.MaxArea;
%eccentricity represents an object's deviation from being a circle
MaxEccentricity=.9999; %0.99
MinEccentricity=PARAMS.MinEccentricity;
%how filled in is the object?
MinSolidity=.6; %0.7
%how bright is the brightest pixel?
minMax=PARAMS.MinMax;
%calculate background by taking an average over 100x100 square
b=imfilter(poi_channel,ones(100)/10000,'symmetric'); %used later to reduce noise from the Smo staining to eliminate pixel values which are unrepresentative of their surroundings.

%calculate properties for each object
possb_cilia_info=regionprops(L,'Eccentricity','Area','Solidity','PixelIdxList', 'BoundingBox', 'MajorAxisLength'); %calculate stats on each 'connected object' identified

%Go through the objects and determine if each is a cilium based on its
%properties.
true_cilia_postn=L; %matrix that tracks the positions of the true cilia
for i=1:length(possb_cilia_info)
    possb_cilia_info(i).positiveCilia = 1; %initialize the positiveCilia feature for the object. This meant to track if the object is a true cilium 
    maxI=max(cilia_channel(possb_cilia_info(i).PixelIdxList)); %get pixel with max intensity from cilia staining for this possible "cilia"
    possb_cilia_info(i).maxI = maxI;
    %first check if any other objects are nearby to weed out any
    %"cilium-like" objects created from noisy acetylated-tubulin staining.
    boundingBox = possb_cilia_info(i).BoundingBox;
    %calculate the x coordinate for the upper left corner of expanded
    %bounding box
    upperLeftX = boundingBox(1) - PARAMS.perimeterThres;
    if upperLeftX < 1
        upperLeftX = 1;
    end
    %calculate the y coordinate for the upper left corner of expanded bounding
    %box
    upperLeftY = boundingBox(2) - PARAMS.perimeterThres;
    if upperLeftY < 1
        upperLeftY = 1;
    end
    %calculate the new width of the bounding box
    matSize = size(true_cilia_postn);
    newWidth = boundingBox(3) + (2*PARAMS.perimeterThres);
    if (upperLeftX + newWidth) > matSize(1)
        newWidth = matSize(1) - upperLeftX;
    end
    %calculate the new length of the bounding box
    newLength = boundingBox(4)+ (2*PARAMS.perimeterThres);
    if (upperLeftY + newLength) > matSize(2)
        newLength = matSize(2) - upperLeftY;
    end
    %determine the values needed for the expanded bounding box
    exp_startx = floor(upperLeftX);
    exp_starty = floor(upperLeftY);
    exp_endx = floor(upperLeftX+newWidth);
    exp_endy = floor(upperLeftY+newLength);
    expandedBox = true_cilia_postn(exp_starty:exp_endy, exp_startx:exp_endx);
    elements = numel(unique(expandedBox)); %get the number of unique values in this expanded box
    if elements > 2 %if there is more than 2 unique values in the expanded box, throw out this object because that means there are 2 or more objects in this box
        true_cilia_postn(possb_cilia_info(i).PixelIdxList)=0; %if not, change this "cilia's" mapping number to 0
        possb_cilia_info(i).positiveCilia = 0; %set cilia marker to 0/FALSE
    %check if this "cilium" posseses the feature/meets the thresholds for
    %being labeled a "cilium"
    elseif maxI<minMax | possb_cilia_info(i).Area<MinArea | possb_cilia_info(i).Area>MaxArea | possb_cilia_info(i).Eccentricity>MaxEccentricity | possb_cilia_info(i).Eccentricity<MinEccentricity | possb_cilia_info(i).Solidity<MinSolidity
        true_cilia_postn(possb_cilia_info(i).PixelIdxList)=0; %if not, change this "cilia's" mapping number to 0
        possb_cilia_info(i).positiveCilia = 0; %set cilia marker to 0/FALSE
    end
end

display('done finding true cilia');

true_cilia_indices=[]; %vector to hold the indices of the true "cilia" as identified by the "positiveCilia" variable
%for each true cilia, try to find Gli2 signal in the Gli2 channel
target_channel_size = size(poi_channel);
%determine the number of true cilia to quantify
ciliaLogical = [possb_cilia_info.positiveCilia]; %logical vector indicating which objects are cilia
allCilia = sum(ciliaLogical); %by summing the logical vector, we calculate the number of true cilia
positiveCiliaCounter = 0; %counter variable to keep track what number cilium were on
poi_intensities = []; %vector to hold the Gli2 intensities from the cilia
poi_areas = []; %vector to hold the area quantified for the POI in each cilium
%calculate background in the Gli2 channel by taking an average over 100x100 square
b=imfilter(poi_channel,ones(100)/10000,'symmetric'); 
%Loop through all possible objects. If the object was identified as a
%cilium, determine the correct POI signal, calculate the mean fluorescence
%from the Gli2 channel, and subtract the background.
for i=1:numel(possb_cilia_info)
    if possb_cilia_info(i).positiveCilia == 1
        true_cilia_indices(end+1) = i;
        positiveCiliaCounter = positiveCiliaCounter + 1;
        boundingBox = possb_cilia_info(i).BoundingBox; %get the bounding box for the cilia
        %get the coordinates of the bounding box for this cilium to create
        %the snapshot picture of the identified cilium
        startx = floor(boundingBox(1));
        starty = floor(boundingBox(2));
        endx = startx + boundingBox(3);
        endy = starty + boundingBox(4);
        %a couple of checks to make sure the snapshot pictures of each
        %cilium print out ok in the tif file
        if startx == 0
            startx = 1;
        end
        if starty == 0
            starty = 1;
        end
        %print out the snapshot in a tif file
        ciliaWindow(positiveCiliaCounter).box = (cilia_channel(starty:endy, startx:endx));
        fig1=figure(1);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter);
        imshow(true_cilia_postn(starty:endy, startx:endx));

        %expand the bounding box to see if we can get Gli2 staining in box
        %calculate the x coordinate for the upper left corner of expanded bounding box
        upperLeftX = boundingBox(1) - (0.3*boundingBox(3));
        if upperLeftX <= 1
            upperLeftX = 1;
        end
        %calculate the y coordinate for the upper left corner of expanded bounding
        %box
        upperLeftY = boundingBox(2) - (0.3*boundingBox(4));
        if upperLeftY <= 1
            upperLeftY = 1;
        end
        %calculate the new width of the bounding box
        newWidth = 1.60*boundingBox(3);
        if (upperLeftX + newWidth) > target_channel_size(1)
            newWidth = target_channel_size(1) - upperLeftX;
        end
        %calculate the new length of the bounding box
        newLength = 1.60*boundingBox(4);
        if (upperLeftY + newLength) > target_channel_size(2)
            newLength = target_channel_size(2) - upperLeftY;
        end
        %determine the coordinate of the expanded bounding box in Gli2
        %channel
        exp_box_startx = floor(upperLeftX);
        exp_box_starty = floor(upperLeftY);
        exp_box_endx = floor(upperLeftX+newWidth);
        exp_box_endy = floor(upperLeftY+newLength);
        %save a snapshot of the expanded bounding box in the Gli2 channel in a tif file
        fig2=figure(2);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter);
        imshow(poi_channel(exp_box_starty:exp_box_endy, exp_box_startx:exp_box_endx));

        %find the right Gli2 signal in the expanded bounding box and
        %quantify it
        poiBox = poi_channel(exp_box_starty:exp_box_endy, exp_box_startx:exp_box_endx);
        %identify pixels with real signal in the Gli2 channel using the
        %user-define threshold
        poiSignal=(medfilt2(poiBox,[3 3])>PARAMS.TargetProtThreshold);
        %show the mask in the Gli2 channel for the expanded
        %bounding box and save a snapshot in a tif file
        fig3=figure(3);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter);
        imshow(poiSignal);
        %determine the objects in Gli2 channel
        targetObj=bwlabel(poiSignal);
        %get info on each possible object in the expanded bounding box in
        %the POI channel
        possb_target_info=regionprops(targetObj,'Eccentricity','Area','Solidity','PixelIdxList', 'BoundingBox', 'MajorAxisLength');
        %loop through objects and determine the one with the greatest
        %intensity and maximum area
        maxPixelIndex = 1;
        maxArea = -1;
        maxPixel = -1;
        for j=1:numel(possb_target_info)
            maxI=max(poiBox(possb_target_info(j).PixelIdxList));
            if maxI > maxPixel 
                if possb_target_info(j).Area > maxArea
                    maxPixelIndex=j;
                    maxPixel = maxI;
                    maxArea = possb_target_info(j).Area;
                end
            end
        end
        %check if we found an object that we can say is the Gli2 signal in
        %this cilium
        if maxPixel > -1 %need this check because sometimes no objects are found in the expanded bounding box because the Gli2 signal is too weak
            %get matrix for background correction
            targetB = b(exp_box_starty:exp_box_endy, exp_box_startx:exp_box_endx);
            %calculate the mean intensity for the Gli2, subtract
            %the background, and store the intensity value
            poi_intensities(positiveCiliaCounter) = mean(poiBox(possb_target_info(maxPixelIndex).PixelIdxList)) - mean(targetB(possb_target_info(maxPixelIndex).PixelIdxList));
            %store the area quantified for the Gli2
            poi_areas(positiveCiliaCounter) = possb_target_info(maxPixelIndex).Area;
            %get the coordinates of the Gli2 quantified region so a snapshot can
            %be saved
            target_startx = floor(possb_target_info(maxPixelIndex).BoundingBox(1));
            target_starty = floor(possb_target_info(maxPixelIndex).BoundingBox(2));
            target_endx = target_startx + possb_target_info(maxPixelIndex).BoundingBox(3);
            target_endy = target_starty + possb_target_info(maxPixelIndex).BoundingBox(4);
            %a couple of checks to make sure pictures print out ok
            if target_startx == 0
                target_startx = 1;
            end
            if target_starty == 0
                target_starty = 1;
            end
        else %if no Gli2/POI signal found
            poi_intensities(positiveCiliaCounter) = 0;
            poi_areas(positiveCiliaCounter) = 0;
            targetSize = size(poiSignal);
            target_startx = 1;
            target_starty = 1;
            target_endx = targetSize(2)-1;
            target_endy = targetSize(1)-1;
        end
        %save a snapshot of the quantified Gli2/POI region in the cilium
        fig4=figure(4);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter);
        imshow(poiSignal(target_starty:target_endy, target_startx:target_endx));
    end
end

saveas(fig1, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_chosen.tif']);
saveas(fig2, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_channel.tif']);
saveas(fig3, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_mask.tif']);
saveas(fig4, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_chosen.tif']);
close(fig1);
close(fig2);
close(fig3);
close(fig4);


end