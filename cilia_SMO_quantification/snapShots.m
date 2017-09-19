%This function outputs two tif files: one contains snapshots for each of
%the true-cilia determined in the "quant_cilia_for_field.m" function and
%the other contains snapshots of what the signal from the POI channel
%looks like in each of the true-cilia. These files are for the user to
%verify that the pipeline is really identifying cilia and quantifying the
%correct POI signal.
function ciliaWindow = snapShots(R, L1,J1,V1, PARAMS)
ciliaWindow = [];
numShots = numel(R); %get number of possible cilium-objects
ciliaLogical = [R.positiveCilia];
% out of R, find positive cilia 
allCilia = sum(ciliaLogical); %get the number of true-cilia
if allCilia == 0 %check to prevent error if no cilia are found
    disp(['Zero cilia in field_' num2str(PARAMS.currentField, '%d')]);
    saveas(figure(1), [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_chosen.tif']);
    saveas(figure(2), [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_.tif']);
    return;
end

positiveCiliaCounter = 0;
%initialize counter
ciliaWindow(allCilia).box = [];
for M = 1 : numShots
    if R(M).positiveCilia == 1 
        positiveCiliaCounter = positiveCiliaCounter + 1; %add to counter
        %get the coordinates of the snapshot
        startx = floor(R(M).BoundingBox(1));
        starty = floor(R(M).BoundingBox(2));
        endx = startx + R(M).BoundingBox(3);
        endy = starty + R(M).BoundingBox(4);
        %a couple of checks to make sure pictures print out ok
        if startx == 0
            startx = 1;
        end
        if starty == 0
            starty = 1;
        end
        ciliaWindow(positiveCiliaCounter).box = (J1(starty:endy, startx:endx));
        fig1=figure(1);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter)
        imshow(L1(starty:endy, startx:endx))
        fig2=figure(2);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter)
        imshow(J1(starty:endy, startx:endx));
    end
end

%save the figures
saveas(fig1, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_chosen.tif']);
saveas(fig2, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_.tif']);
close(fig1);
close(fig2);

end