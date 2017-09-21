%specify the path to the directory holding the scripts from the "cilia_SMO_quantification" folder
cd('replace_with_path_to_the_scripts_from cilia_SMO_quantification_folder');

%PARAMETERS that can be altered based on the imaging
PARAMS.Threshold=95;%This is the threshold for considering a pixel as signal or noise in the cilia staining channel
PARAMS.MinMax=190; %This is the threshold for the minimum maximum intensity value for a potential "cilium" object to be considered as a true cilium
PARAMS.MinArea=21;%This is the threshold for the minimum area for a potential-cilium" object to be considered as a true cilium. This may have to change based on the cilium-shape and resolution of the image.
PARAMS.MaxArea=200;%This is the threshold for the maximum area for a potential-cilium" object to be considered as a true cilium. This may have to change based on the cilium-shape and resolution of the image.
PARAMS.MinEccentricity=.85; %This is the threshold for the minimum eccentricity for a cilium object. This threshold may change depending on the cilium-shape in different cell lines/types.
%Specify the order of the images in the Leica (.lif) file
PARAMS.dapiNum=1; %to specify the position of the DAPI channel
PARAMS.ciliaNum=2; %to specify the position of the cilia-marker channel
PARAMS.targetNum=3; %to specify the position of the channel for the POI
PARAMS.perimeterThres = 5; % Distance to search around each potential-cilium object for other objects

PARAMS.numFields=3; %to specify the number of fields captured in this Leica file
PARAMS.numMarkers=3; %to specify the number of markers searched for; usually just 3 because one staining each for DAPI, cilia, POI
%specify the leica file to analyze
leicaFile = 'replace_with_name_of_Leica_file';
%generate output file name
PARAMS.output_name=leicaFile(1:end-4);
%get the max-Z projections for each marker
[fields] = get_max_z_projections(leicaFile, PARAMS);

%Begin quantifying cilia and POI in cilia
true_cilia_area = []; %vector to hold the area of each of cilia
true_cilia_length = []; %vector to hold the length of each of cilia
true_poi_intensities = []; %vector to hold the signal for the POI in each cilium
for(i = 1:PARAMS.numFields) %go through each field to quantify POI in each cilium
    PARAMS.currentField=i;
    %call to function to actually quantify
    [true_cilia_postn,cilia_channel,poi_channel,poi_intensities,true_cilia_indices,possb_cilia_info]=quant_cilia_for_field(fields{i,1}, PARAMS);
    %call to function to output snapshots of each chosen cilium and the POI
    %signal in each cilium
    ciliaWindow = snapShots(possb_cilia_info, true_cilia_postn,poi_channel,poi_intensities, PARAMS);
    
    %aggregate data from different fields
    for i = 1:numel(true_cilia_indices)
        index = true_cilia_indices(i);
        true_poi_intensities(end+1) = poi_intensities(i);
        true_cilia_area(end+1) = possb_cilia_info(index).Area;
        true_cilia_length(end+1) = possb_cilia_info(index).MajorAxisLength;
    end
end
disp(['Quantified ',num2str(length(true_poi_intensities)), ' cilia']);
%output statistics on each identified cilia
vals = [true_poi_intensities; true_cilia_area; true_cilia_length]';
csvwrite([PARAMS.output_name '.csv'], vals);