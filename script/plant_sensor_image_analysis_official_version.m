%%
clear all
close all

% USER IMPORTANT INFORMATION: 5 stars ***** means that the code needs adjustment (correct path, filenames, etc) !!!!!!!

% ***** run this first, paths and file names of all file of interest for this script
% path to the chlorophyll files
pathpath = 'C:\projekti\projekti NIB\lezije\CompleteAnalysis\Experiment_2\L2\ch00_chlorophyll\';
Files = dir(strcat(pathpath,'*.tif'));
Files_405 = dir('C:\projekti\projekti NIB\lezije\CompleteAnalysis\Experiment_2\L2\ch01_405\*.tif');
Files_488 = dir('C:\projekti\projekti NIB\lezije\CompleteAnalysis\Experiment_2\L2\ch02_488\*.tif');
% load the excel file and read in the positions of each lesion for each ROI01
% the files in the excel should be named the same as the files in the folders, but without the ending Processed_ch00 etc 
positions = readtable('C:\projekti\projekti NIB\lezije\CompleteAnalysis\Experiment_2\SpatialAnalysis\RywalNahG_2nd_exp_lesion_position_L2.xlsx');

%% 
% plant_sensor_image_analysis_official_version.m
%
% What does the script do?
% The scripts uploads a series of images from the confocal microscope,
% taken at three different wavelengths, chlorophyll, 405nm and 488 nm. It
% then segments out individual chloroplasts and calculates the 405/488
% ratio in the chloroplasts  - this is  measure of the redox state in the
% tissue. It also imports information about the position of the lesion for
% each image and calculates the spatial distribution of the redox state. 
%
% Date of last change: November 22, 2020
% Update: clean code + automatic recognition of image size
% Author: Anze Zupanic

%% DATA AND IMAGE PREPARATION FOR THE ANALYSIS
% to perform all the different steps of the analysis in this script, the
% following is required:
    % microscopic images, 3 different channels, chlorophyl channel, 405 nm
    % channel and 488 nm channel. IMPORTANT: there should be a separate
    % folder for each type of images, and the names of the files should be
    % the SAME, except for the ending of the files 
    % These are good names: 
        % 17022020_Rywal_NahG_pt-roGFP_L2_P1_L2_PVY-WILGA_3dpi.lif_lesion01_ROI01_Processed001_ch00.tif
        % 17022020_Rywal_NahG_pt-roGFP_L2_P1_L2_PVY-WILGA_3dpi.lif_lesion01_ROI01_Processed001_ch01.tif
        % 17022020_Rywal_NahG_pt-roGFP_L2_P1_L2_PVY-WILGA_3dpi.lif_lesion01_ROI01_Processed001_ch02.tif
    % For the spatial analysis, each file name (only one channel can be
    % used here, it's just important that the alphabetical order matches),
    % should be accompanied by the position of the lesion
% Example files, with appropriate naming and formats can be found on 
% https://github.com/NIB-SI/SensorPlantAnalysis

%% PART 1: BATCH IMPORT OF IMAGES AND DETECTION OF CHLOROPLASTS
    % extracts the chloroplast and creates a mask
    % Outputs 
    %   Images: image diary depicting sequential image operations made
    %   Variables:
        % masks - a binary image with signal only present within detect
            % chloroplasts
        % chloroplast_pixels - which pixels belong to each chloroplast in each image 
        % chloroplastStats - centroid coordinates, area and
            % mean_fluorescence intensity for each chloroplast in each
            % image
        % figureStats - for each image, number of chloroplasts, their
            % average intensity and total area of chloroplasts
        % Files - names of the files analyzed         
        
% reads all tiff files
Files = dir(strcat(pathpath,'*.tif'));
numFiles = numel(Files); % returns number of files found

% loop over all files
for jj = 1:numFiles
    jj 
    % read in each image and convert to grayscale
    image_original = imread(strcat(Files(jj).folder,'\',Files(jj).name));
    image_gray = rgb2gray(image_original);
    
    % first a very simplistic way of getting rid of noise (background signal
    % from other cells,etc). Finds the peak of the histogram and deletes
    % everything to the right of the peak +some offset
    [counts,x] = imhist(image_gray);
    [i,whereMax] = max(counts);
    % ***** can try other values here
    image_gray(image_gray<whereMax+5)=0;
    
    % make a binary image, based on a local adaptive method
    % parameters:
    %   sensitivity = 0.6 (higher sensitivity means more noise gets in,
        %   lower means less signal, you need to find a good compromise). This
        %   needs to be tested for each batch of images. 
    %   neighbourhood size = 11 (lower values separate better between chloroplasts)
    % ***** (try playing with settings a bit for each new microscope dataset)
    T = adaptthresh(image_gray,0.6,'NeighborhoodSize',11, 'Statistic', 'Gaussian');
    image_binary = imbinarize(image_gray, T);
    
    % erosion, dilation, median filtering and areaopening operations provide
    % good results, regardless of image quality
    se = strel('square', 6); % ***** other shapes and sizes can be tried for the morphological structure element
    image_open = imopen(image_binary,se);
    image_filtered = medfilt2(image_open);
    % here the minimum size of the chloroplast is defined, anything smaller
    % gets deleted
    % *****define the minimum size according to species/imaging protocol
    image_sized = bwareaopen(image_filtered, 60);
    % OUTPUT: final mask for chloroplasts in this image:
    masks{jj} = image_sized;
    
    % calculate total chloroplast are in the image (area of mask)
    chloroplastArea{jj} = bwarea(image_sized);

    % NORMALIZATION
    % calculate mean fluorescence of chloroplasts
    %   - change image to double * 255 gives the same values as in image_gray
    %   - then multiply with mask, now you have 0s where no mask and gray 
    %       values where there is a mask
    %   - mean(mean*image_temp)) now gives the average value of the whole
    %       image (and not only mask), so it needs to be normalized 
    image_temp = im2double(image_gray).*255;
    image_temp = image_temp.*image_sized;
    meanChloroplastFluorescence{jj} = numel(image_temp(:))*mean(mean(image_temp))/chloroplastArea{jj};
    
    % get the individual chloroplasts or clusters of chloroplasts 
    % where it-s difficult to divide them into individual ones.
    % chloroplast_pixels.PixelldxList - indexes where the pixels belonging to
    %   chloroplasts are (careful, matlab indexes by rows!!!)
    % chloroplast_pixels{jj}.NumObjects / number of chloroplasts identified
    % OUTPUT - pixel indexes belonging to individual chloroplasts
    chloroplast_pixels{jj} = bwconncomp(image_sized);
    numChloroplastElements(jj) = chloroplast_pixels{jj}.NumObjects;

    % prepare for visualization of individual chloroplasts
    labeledChloroplasts = labelmatrix(chloroplast_pixels{jj});
    RGB_label = label2rgb(labeledChloroplasts,'spring','c','shuffle');

    % regionprops measures properties for each connected object in output
    % of bwconncomp. Here, it-s the Area and centroid for each
    % chloroplasts and the MeanIntensity 
    grayscaleStats = regionprops(chloroplast_pixels{jj},image_gray,{'Area','Centroid', 'MeanIntensity'});
    % OUTPUT here the stats per individual chloroplast are saved
    chloroplastStats{jj} = grayscaleStats;
    
    % put all figure stats into a single table
    figureStats(jj,1) = jj;
    figureStats(jj,2) = numChloroplastElements(jj);
    figureStats(jj,3) = chloroplastArea{jj};
    figureStats(jj,4) = meanChloroplastFluorescence{jj};

    % PLOTTING
    % a single figure that illustrates the image analysis process
    % ***** uncomment this if you want the figures of the individual steps of chloroplast recognition saved 
    FigH = figure('Position', get(0, 'Screensize'));
    set(FigH, 'Visible', 'off');
    subplot(3,2,1)
    imshow(image_original)
    title('Original image')
    subplot(3,2,2)
    imshow(image_gray)
    % ***** the top value here depends on the intensity of the images,
    % choose a value high enough (close to max intensity of all images)
    caxis([0 40])
    title('Grayscale Image Scaled')
    subplot(3,2,3)
    imshow(image_binary)
    title('After binarization')
    subplot(3,2,4)
    imshow(image_sized)
    title('After filtering and morphological operations')
    subplot(3,2,5)
    imshowpair(image_binary,image_sized)
    title('Before after filtering difference')
    subplot(3,2,6)
    imshow(RGB_label)
    title('Mean fluorescence of each connected chloroplast');
    hold on
    for kk = 1 : numChloroplastElements(jj)
        text(grayscaleStats(kk).Centroid(1),grayscaleStats(kk).Centroid(2), ...
            sprintf('%2.1f', grayscaleStats(kk).MeanIntensity), ...
            'EdgeColor','b','Color','k');
    end
    hold off
    % to show the figure
    % figure(FigH)
    
    % select path and save
    % we select the printing resolution
    iResolution = 300;
    % we select to crop or not the figure
    set(gcf, 'PaperPositionMode', 'auto');
    % saving full screen size
    F    = getframe(FigH);
    imwrite(F.cdata, strcat(pathpath, Files(jj).name, '.jpeg'), 'jpeg')
end

% OUTPUT here are all the stats of the individual figures
figureStats = array2table(figureStats, 'VariableNames', {'figureNumber', 'NumberChloroplasts', 'AreaChloroplasts','meanChloroplastFluorescence'});
% saving the important results
savefile = strcat(pathpath,'analyzedChloroplasts.mat');
% here we define all OUTPUTs to save
save(savefile, 'masks','chloroplast_pixels','chloroplastStats', 'figureStats', 'Files');

%% load results of segmentation
load(strcat(pathpath,'analyzedChloroplasts.mat'))


%% PART 2: 405/488 RATIO per FIGURE and per CHLOROPLAST
    % masks from first part are used on the 405 and 488 microscopic images
    % the ratio between the two fluorescences is calculated
        
numFiles = numel(Files_405);
% read the images one by one and do the analysis
for jj = 1:numFiles
    jj
    % read in the image and convert to grayscale for the 405 nm channels
    strcat(Files_405(jj).folder,'\',Files_405(jj).name);
    image_405 = imread(strcat(Files_405(jj).folder,'\',Files_405(jj).name));
    image_405_gray = rgb2gray(image_405);
    % convert to double so that multiplying with mask is possible
    image_405_double = im2double(image_405_gray).*255;
    
    % same for 488
    strcat(Files_488(jj).folder,'\',Files_488(jj).name);
    image_488 = imread(strcat(Files_488(jj).folder,'\',Files_488(jj).name));
    image_488_gray = rgb2gray(image_488);
    image_488_double = im2double(image_488_gray).*255;
    
    % calculate the ratio between the images, without any changes 405/488
    % division comes first, then the application of mask
    image_rat = image_405_double./image_488_double;
    image_ratio = image_rat.*masks{jj};
    % OPTIONAL image: to show what the ratio looks like
%     figure
%     imshowpair(image_1, image_12,'montage')
%     colorbar
    
    % OUTPUT all the image ratios are saved
    ratios{jj} = image_ratio;
    
    % NORMALIZATION/OUTPUT
    % here the mean ratios per figure are saved
    % calculate the mean ratio, again we need to normalize only to the part
    % of the image that is the chloroplast (explanation above in PART 1)
    normalizedRatio{jj} = numel(image_ratio(:))*mean(mean(image_ratio))/figureStats.AreaChloroplasts(jj);
    
    
    % PLOTTING
    % series of images to illustrate, what was done
    % ***** uncomment this if you want plotting and saving the images
    FigH = figure('Position', get(0, 'Screensize'));
    set(FigH, 'Visible', 'off');
    subplot(2,2,1)
    imshow(image_405_gray)
    caxis([0 30])
    title('Scaled 405')
    subplot(2,2,2)
    imshow(image_488_gray)
    caxis([0 30])
    title('Scaled 488')
    subplot(2,2,3)
    image(image_rat,'CDataMapping','scaled')
    axis equal
    colormap(gray)
    % ***** check the histogram and then define the color axis
    % hist(image_rat(:),100)
    caxis([0 4])
    title('405/488')
    subplot(2,2,4)
    image(image_ratio,'CDataMapping','scaled')
    axis equal
    colormap(gray)
    % ***** check the histogram and then define the color axis
    % hist(image_ratio(:),100)
    caxis([0 4])
    title('405/488 in chloroplasts')
    % to show the figure
    % figure(FigH)
    
    % select path and save
    % we select the printing resolution
    iResolution = 300;
    % we select to crop or not the figure
    set(gcf, 'PaperPositionMode', 'auto');
    % saving full screen size
    F    = getframe(FigH);
    imwrite(F.cdata, strcat(pathpath, Files(jj).name, '_ratio','.jpeg'), 'jpeg')
     

    % PER CHLOROPLAST
    % OUTPUT
    % her the ratio per all chloroplasts per each figure are saved
    grayscaleChloroplasts = regionprops(masks{jj},image_ratio,{'Centroid', 'MeanIntensity'});
    chloroplast_ratios{jj} = grayscaleChloroplasts;

end

% PLOTTING the results PER FIGURE
% it is very important that the figures that are analyzed are in a sensible
% order, because the points on the plot are in the same order as the
% figures that were analyzed! 
normRatio = cell2mat(normalizedRatio);
figure
plot(normRatio,'o', 'MarkerSize', 20,'MarkerFaceColor', 'b')
set(gca,'FontSize',20)
% ***** here you can put separators between the samples ***** 
% for example, if the first 10 are one treatment and the rest are another
% treatment, you can put lines in the correct place
% hold on
% xline(10.5)
% OR you can color the data points according to treatments
% plot(1:10,normRatio([1:5,6:10]),'or','MarkerSize',12, 'MarkerFaceColor', 'r')
%   hold on
% plot(11:20,normRatio([11:20]),'ob','MarkerSize',12, 'MarkerFaceColor', 'b')
%   the first then would be red color in this case


% PLOTTING the results PER CHLOROPLAST   
% plot without ordering, plot a distribution over all chloroplasts as
% boxplots
% ***** important, put the bplot.m script in the directory in which you are
figure
hold on
for jj = 1:numFiles
   temp = [chloroplast_ratios{jj}.MeanIntensity];
   bplot(temp,jj)
end
set(gca,'FontSize',20)
hold off

% saving the important results of the 405/488 analysis
savefile = strcat(pathpath,'chloroplastRatios.mat');
% here we define all OUTPUTs to save
save(savefile, 'chloroplast_ratios', 'normalizedRatio');


%% load results of 405/488 analysis
load(strcat(pathpath,'chloroplastRatios.mat'))


%% PART 3A: SPATIAL ANALYSIS INCLUDED FOR THE REGIONS OF INTEREST CLOSEST TO THE LESION (ROI1)

numFiles = numel(Files_405);
% find all images, which have ROI1 in their name and not MOCK, DTT, h202 or
% control
index = [];
for jj = 1:size(Files,1) 
   if contains(Files_405(jj).name,'ROI01') && ~contains(Files_405(jj).name,{'MOCK','DTT','H2O2','control', 'ddH2O'})
   % here for ROI2
 % if contains(Files_405(jj).name,'ROI02') && ~contains(Files_405(jj).name,{'MOCK','DTT','H2O2','control', 'ddH2O'})
     index = [index,jj];             
   end  
end
Files_filtered = Files(index,:);
Files_405_filtered = Files_405(index,:);
Files_488_filtered = Files_488(index,:);
chloroplastStats_filtered = chloroplastStats(1,index);
figureStats_filtered = figureStats(index,:);
masks_filtered = masks(index);


% we will now add the distance to the lesion to the chloroplastStats
% variable
% ***** size of the image in this case is 512 pixels, for images of different
% size, write a different number
% for each image
% we repeat the same analysis as before, but with spatial information
for jj = 1:size(Files_filtered,1)
    % in this first part the intensities per chloroplast are calculated as
    % before
    
    % PER FIGURE
    % read in the image and convert to grayscale
    strcat(Files_405_filtered(jj).folder,'\',Files_405_filtered(jj).name);
    image_405 = imread(strcat(Files_405_filtered(jj).folder,'\',Files_405_filtered(jj).name));
    image_405_gray = rgb2gray(image_405);
    % convert to double so that multiplying with mask is possible
    image_405_double = im2double(image_405_gray).*255;
    
    % same for 488
    strcat(Files_488_filtered(jj).folder,'\',Files_488_filtered(jj).name);
    image_488 = imread(strcat(Files_488_filtered(jj).folder,'\',Files_488_filtered(jj).name));
    image_488_gray = rgb2gray(image_488);
    image_488_double = im2double(image_488_gray).*255;
    
    % calculate the ratio between the images, without any changes 405/488
    % division comes first, then the application of mask
    image_rat = image_405_double./image_488_double;
    image_ratio = image_rat.*masks_filtered{jj};
    % if removal of background uncomment the next two lines and comment the
    % previous two
%     image_rat = image_405_double./image_488_removedBackground;
%     image_ratio = image_rat.*masks_filtered{jj};
    % this is to show what the ratio looks like
%     figure
%     imshowpair(image_1, image_12,'montage')
%     colorbar
    
    % OUTPUT all the image ratios are saved
    ratios{jj} = image_ratio;
    
    % NORMALIZATION/OUTPUT
    % here the mean ratios per figure are saved
    % calculate the mean ratio, again we need to normalize only to the part
    % of the image that is the chloroplast (explanation above)
    normalizedRatio_background{jj} = numel(image_ratio(:))*mean(mean(image_ratio))/figureStats_filtered.AreaChloroplasts(jj);
    
    % PER CHLOROPLAST
    % OUTPUT
    % her the ratio per all chloroplasts per each figure are saved
    grayscaleChloroplasts = regionprops(masks_filtered{jj},image_ratio,{'Centroid', 'MeanIntensity'});
    chloroplast_ratios_background{jj} = grayscaleChloroplasts;
    
    % in this second part, the distance between each chloroplast and the
    % lesion is calculated
    
    numChloroplasts = size(chloroplastStats_filtered{1,jj});
    % for each chloroplast
    for kk = 1:numChloroplasts(1)
        % calculate the distance to the lesion
        coordinatesChloroplast = chloroplastStats_filtered{1,jj}(kk).Centroid;
        switch positions{jj,2}{1}
            case 'spodaj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = 512 - coordinatesChloroplast(2);
            case 'zgoraj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = coordinatesChloroplast(2);
            case 'levo'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = coordinatesChloroplast(1);
            case 'desno'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = 512 - coordinatesChloroplast(1);
            case 'levo spodaj'
                % if stat toolbox available change to pdist([coordinatesChloroplast;0,0])
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((coordinatesChloroplast(1))^2+(512 - coordinatesChloroplast(2))^2);
            case 'desno spodaj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((512 -coordinatesChloroplast(1))^2+(512 - coordinatesChloroplast(2))^2);
            case 'levo zgoraj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((coordinatesChloroplast(1))^2+(coordinatesChloroplast(2))^2);
            case 'desno zgoraj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((512 -coordinatesChloroplast(1))^2+(coordinatesChloroplast(2))^2);
        end
        % ***** assign a bin to each chloroplast, bin 1 is closest to
        % lesion, while bin 5 is furthest. If you want more bins, change
        % the code accordingly. 
        if chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 400
            chloroplast_ratios_background{1,jj}(kk).Bin = 5;
        elseif chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 300 && chloroplast_ratios_background{1,jj}(kk).DistanceLesion < 400
            chloroplast_ratios_background{1,jj}(kk).Bin = 4;
        elseif chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 200 && chloroplast_ratios_background{1,jj}(kk).DistanceLesion < 300
            chloroplast_ratios_background{1,jj}(kk).Bin = 3;
        elseif chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 100 && chloroplast_ratios_background{1,jj}(kk).DistanceLesion < 200
            chloroplast_ratios_background{1,jj}(kk).Bin = 2;
        else
            chloroplast_ratios_background{1,jj}(kk).Bin = 1;
        end
    end
end


% PLOT results for each image   
% boxplots
% ***** important, put the bplot.m script in the directory in which you are
IntensityMeans_ROI1 = zeros(size(Files_filtered,1),5);
IntensityVariability_ROI1 = zeros(size(Files_filtered,1),5);
for jj = 1:size(Files_filtered,1)
    jj
    tempStruct = chloroplast_ratios_background{jj};
    % this is a very important line of code that shows how to extract a
    % field from s structure (inside parentheses)
    % in temp vector, for each chloroplast in the image, it is written in
    % which bin it belongs
    tempVector = [tempStruct.Bin];
    % fluorescence per bin per image
    FigH = figure('Position', get(0, 'Screensize'));
    set(FigH, 'Visible', 'off');
    hold on
    for kk = 1:5 % **** because five bins, change for other numbers
       index = tempVector == kk; 
       temp = [tempStruct.MeanIntensity];
       temp = temp(index);
       IntensityMeans_ROI1(jj,kk)=mean(temp);
       IntensityVariability_ROI1(jj,kk)=var(temp);
       bplot(temp,kk); %again 5 bins *****
    end
    set(gca,'FontSize',20)
    hold off
    iResolution = 300;
    % we select to crop or not the figure
    set(gcf, 'PaperPositionMode', 'auto');
    % saving full screen size
    F    = getframe(FigH);
    imwrite(F.cdata, strcat(pathpath, Files_filtered(jj).name, '_spatial_ROI1', '.jpeg'), 'jpeg')
end

% saving the important results
savefile = strcat(pathpath,'spatialAnalysis_ROI1.mat');
% here we define all OUTPUTs to save
save(savefile, 'IntensityMeans_ROI1', 'IntensityVariability_ROI1');

%% PART 3B: SPATIAL ANALYSIS INCLUDED FOR THE REGIONS OF INTEREST FURTHER FROM THE LESION (ROI2)

numFiles = numel(Files_405);
% find all images, which have ROI2 in their name and not MOCK, DTT, h202 or
% control
index = [];
for jj = 1:size(Files,1) 
   if contains(Files_405(jj).name,'ROI02') && ~contains(Files_405(jj).name,{'MOCK','DTT','H2O2','control', 'ddH2O'})
     index = [index,jj];             
   end  
end
Files_filtered = Files(index,:);
Files_405_filtered = Files_405(index,:);
Files_488_filtered = Files_488(index,:);
chloroplastStats_filtered = chloroplastStats(1,index);
figureStats_filtered = figureStats(index,:);
masks_filtered = masks(index);


% we will now add the distance to the lesion to the chloroplastStats
% variable
% ***** size of the image in this case is 512 pixels, for images of different
% size, write a different number
% for each image
% we repeat the same analysis as before, but with spatial information
for jj = 1:size(Files_filtered,1)
    % in this first part the intensities per chloroplast are calculated as
    % before
    
    % PER FIGURE
    % read in the image and convert to grayscale
    strcat(Files_405_filtered(jj).folder,'\',Files_405_filtered(jj).name);
    image_405 = imread(strcat(Files_405_filtered(jj).folder,'\',Files_405_filtered(jj).name));
    image_405_gray = rgb2gray(image_405);
    % convert to double so that multiplying with mask is possible
    image_405_double = im2double(image_405_gray).*255;
    
    % same for 488
    strcat(Files_488_filtered(jj).folder,'\',Files_488_filtered(jj).name);
    image_488 = imread(strcat(Files_488_filtered(jj).folder,'\',Files_488_filtered(jj).name));
    image_488_gray = rgb2gray(image_488);
    image_488_double = im2double(image_488_gray).*255;
    
    % calculate the ratio between the images, without any changes 405/488
    % division comes first, then the application of mask
    image_rat = image_405_double./image_488_double;
    image_ratio = image_rat.*masks_filtered{jj};
    % if removal of background uncomment the next two lines and comment the
    
    % OUTPUT all the image ratios are saved
    ratios{jj} = image_ratio;
    
    % NORMALIZATION/OUTPUT
    % here the mean ratios per figure are saved
    % calculate the mean ratio, again we need to normalize only to the part
    % of the image that is the chloroplast (explanation above)
    normalizedRatio_background{jj} = numel(image_ratio(:))*mean(mean(image_ratio))/figureStats_filtered.AreaChloroplasts(jj);
    
    % PER CHLOROPLAST
    % OUTPUT
    % her the ratio per all chloroplasts per each figure are saved
    grayscaleChloroplasts = regionprops(masks_filtered{jj},image_ratio,{'Centroid', 'MeanIntensity'});
    chloroplast_ratios_background{jj} = grayscaleChloroplasts;
    
    % in this second part, the distance between each chloroplast and the
    % lesion is calculated
    
    numChloroplasts = size(chloroplastStats_filtered{1,jj});
    % for each chloroplast
    for kk = 1:numChloroplasts(1)
        % calculate the distance to the lesion
        coordinatesChloroplast = chloroplastStats_filtered{1,jj}(kk).Centroid;
        switch positions{jj,2}{1}
            case 'spodaj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = 512 - coordinatesChloroplast(2);
            case 'zgoraj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = coordinatesChloroplast(2);
            case 'levo'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = coordinatesChloroplast(1);
            case 'desno'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = 512 - coordinatesChloroplast(1);
            case 'levo spodaj'
                % if stat toolbox available change to pdist([coordinatesChloroplast;0,0])
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((coordinatesChloroplast(1))^2+(512 - coordinatesChloroplast(2))^2);
            case 'desno spodaj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((512 -coordinatesChloroplast(1))^2+(512 - coordinatesChloroplast(2))^2);
            case 'levo zgoraj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((coordinatesChloroplast(1))^2+(coordinatesChloroplast(2))^2);
            case 'desno zgoraj'
                chloroplast_ratios_background{1,jj}(kk).DistanceLesion = sqrt((512 -coordinatesChloroplast(1))^2+(coordinatesChloroplast(2))^2);
        end
        % ***** assign a bin to each chloroplast, bin 1 is closest to
        % lesion, while bin 5 is furthest. If you want more bins, change
        % the code accordingly. 
        if chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 400
            chloroplast_ratios_background{1,jj}(kk).Bin = 5;
        elseif chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 300 && chloroplast_ratios_background{1,jj}(kk).DistanceLesion < 400
            chloroplast_ratios_background{1,jj}(kk).Bin = 4;
        elseif chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 200 && chloroplast_ratios_background{1,jj}(kk).DistanceLesion < 300
            chloroplast_ratios_background{1,jj}(kk).Bin = 3;
        elseif chloroplast_ratios_background{1,jj}(kk).DistanceLesion >= 100 && chloroplast_ratios_background{1,jj}(kk).DistanceLesion < 200
            chloroplast_ratios_background{1,jj}(kk).Bin = 2;
        else
            chloroplast_ratios_background{1,jj}(kk).Bin = 1;
        end
    end
end


% Same plots for ROI2 
IntensityMeans_ROI2 = zeros(size(Files_filtered,1),5);
IntensityVariability_ROI2 = zeros(size(Files_filtered,1),5);
for jj = 1:size(Files_filtered,1)
    jj
    tempStruct = chloroplast_ratios_background{jj};
    % this is a very important line of code that shows how to extract a
    % field from s structure (inside parentheses)
    % in temp vector, for each chloroplast in the image, it is written in
    % which bin it belongs
    tempVector = [tempStruct.Bin];
    % fluorescence per bin per image
    FigH = figure('Position', get(0, 'Screensize'));
    set(FigH, 'Visible', 'off');
    hold on
    for kk = 1:5 % **** because five bins, change for other numbers
       index = tempVector == kk; 
       temp = [tempStruct.MeanIntensity];
       temp = temp(index);
       IntensityMeans_ROI2(jj,kk)=mean(temp);
       IntensityVariability_ROI2(jj,kk)=var(temp);
       bplot(temp,kk); %again 5 bins *****
    end
    set(gca,'FontSize',20)
    hold off
    iResolution = 300;
    % we select to crop or not the figure
    set(gcf, 'PaperPositionMode', 'auto');
    % saving full screen size
    F    = getframe(FigH);
    imwrite(F.cdata, strcat(pathpath, Files_filtered(jj).name, '_spatial_ROI2', '.jpeg'), 'jpeg')
end

% saving the important results
savefile = strcat(pathpath,'spatialAnalysis_ROI2.mat');
% here we define all OUTPUTs to save
save(savefile, 'IntensityMeans_ROI2', 'IntensityVariability_ROI2', 'Files_filtered', 'Files_405_filtered', 'Files_488_filtered');
%% loading the spatial analysis results
load(strcat(pathpath,'spatialAnalysis_ROI1.mat'))
load(strcat(pathpath,'spatialAnalysis_ROI2.mat'))

%% PART 4 NORMALIZATION AND PLOTTING OF THE WHOLE SPATIAL ANALYSIS 
% Normalization of the data against Bin 5 (most distant from lesion)
% for ROI1
IntensityVariability_Norm_ROI1 = IntensityVariability_ROI1./(IntensityMeans_ROI1(:,5));
IntensityMeans_Norm_ROI1 = IntensityMeans_ROI1./(IntensityMeans_ROI1(:,5));
filename = strcat(pathpath, '_SpatialBins','.xlsx');
writematrix(IntensityMeans_ROI1,filename,'Sheet',2)
writematrix(IntensityMeans_Norm_ROI1,filename,'Sheet',3)
% mean, variability and CI for fluorescence intensities (ratios)
MeanIntensityMeans_ROI1 = mean(IntensityMeans_ROI1,'omitnan');
SEMIntensityMeans_ROI1 = std(IntensityMeans_ROI1,'omitnan')/sqrt(size(Files_filtered,1));
CI95IntensityMeans_ROI1 = [MeanIntensityMeans_ROI1+1.96*SEMIntensityMeans_ROI1;
    MeanIntensityMeans_ROI1-1.96*SEMIntensityMeans_ROI1];

% for ROI2
IntensityVariability_Norm_ROI2 = IntensityVariability_ROI2./(IntensityMeans_ROI2(:,5));
IntensityMeans_Norm_ROI2 = IntensityMeans_ROI2./(IntensityMeans_ROI2(:,5));
writematrix(IntensityMeans_ROI2,filename,'Sheet',4)
writematrix(IntensityMeans_Norm_ROI2,filename,'Sheet',5)
% mean, variability and CI for fluorescence intensities (ratios)
MeanIntensityMeans_ROI2 = mean(IntensityMeans_ROI2,'omitnan');
SEMIntensityMeans_ROI2 = std(IntensityMeans_ROI2,'omitnan')/sqrt(size(Files_filtered,1));
CI95IntensityMeans_ROI2 = [MeanIntensityMeans_ROI2+1.96*SEMIntensityMeans_ROI2;
    MeanIntensityMeans_ROI2-1.96*SEMIntensityMeans_ROI2];




% Here the panel figure, in the first row spatial fluorescence for ROI1 and
% 3dpi, 5dpi, 7dpi. In the bottom row the same for ROI2.
% first ROI1, dpi3
% ########################
FigH = figure('Position', get(0, 'Screensize'));
set(FigH, 'Visible', 'off');
iFontSize = 28;
strFontUnit = 'points'; % [{points} | normalized | inches | centimeters | pixels]
strFontName = 'Arial'; % [Times | Courier | OTHERS ]
strFontWeight = 'normal'; % [light | {normal} | demi | bold]
strFontAngle = 'normal'; % [{normal} | italic | oblique] ps: only for axes
fLineWidth = 4.0; % width of the line of the axes

subplot_er(2,3,1)
index3dpi = [];
for jj = 1:size(Files_filtered,1) 
   if contains(Files_405_filtered(jj).name,'3dpi')
     index3dpi = [index3dpi,jj];             
   end  
end
IntensityMeans_Norm_3dpi_ROI1 = IntensityMeans_Norm_ROI1(index3dpi,:);
MeanIntensityMeans_Norm_3dpi_ROI1 = mean(IntensityMeans_Norm_3dpi_ROI1,'omitnan');
SEMIntensityMeans_Norm_3dpi_ROI1 = std(IntensityMeans_Norm_3dpi_ROI1,'omitnan')/sqrt(size(index3dpi,2));
CI95IntensityMeans_Norm_3dpi_ROI1 = [MeanIntensityMeans_Norm_3dpi_ROI1+1.96*SEMIntensityMeans_Norm_3dpi_ROI1;
    MeanIntensityMeans_Norm_3dpi_ROI1-1.96*SEMIntensityMeans_Norm_3dpi_ROI1];
% first part of figure, the individual points
% first the boxplot
bplot(IntensityMeans_Norm_3dpi_ROI1(:,1),1, 'linewidth', 4)
hold on
bplot(IntensityMeans_Norm_3dpi_ROI1(:,2),2, 'linewidth', 4)
bplot(IntensityMeans_Norm_3dpi_ROI1(:,3),3, 'linewidth', 4)
bplot(IntensityMeans_Norm_3dpi_ROI1(:,4),4, 'linewidth', 4)
% then the individual points
x_axis = repmat([1:4]',size(IntensityMeans_Norm_3dpi_ROI1,1),1);
IntensityMeans_Norm_T_3dpi_ROI1 = IntensityMeans_Norm_3dpi_ROI1(:,1:4)';
y_axis = IntensityMeans_Norm_T_3dpi_ROI1(:);
scatter(x_axis, y_axis, ...
    'MarkerEdgeColor',"#999999",...
              'MarkerFaceColor',"#999999",...
              'LineWidth',4)
line([0 4.75], [1 1],'LineWidth',3,'color',"#999999", 'LineStyle', '--')
% % finally the lines between the points
% plot(x_axis, y_axis, ':',...
%              'Color', "#999999",...
%               'LineWidth',2)
set(gca, ...
    'XTick', 1:1:4, ... ticks of x axis
    'YTick', 0:1:ceil(max(max(IntensityMeans_Norm_ROI1)))-1, ... ticks of y axis
    'XTickLabel', [],...
    'XGrid', 'off', ... [on | {off}]
    'YGrid', 'off', ... [on | {off}]
    'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorGrid', 'off' , ... [on | {off}]
    'YMinorGrid', 'off', ... [on | {off}]
    'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorTick', 'off' , ... [on | {off}]
    'YMinorTick', 'off', ... [on | {off}]
    'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
    'TickLength', [.02 .02], ... length of the ticks
    'XColor', [.1 .1 .1], ... color of x axis
    'YColor', [.1 .1 .1], ... color of y axis
    'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
    'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
    'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'FontName', strFontName, ... kind of fonts of labels
    'FontSize', iFontSize, ... size of fonts of labels
    'FontUnits', strFontUnit, ... units of the size of fonts
    'FontWeight', strFontWeight, ... weight of fonts of labels
    'FontAngle', strFontAngle, ... inclination of fonts of labels
    'LineWidth', fLineWidth, ... % width of the line of the axes
    'xlim',[0.25 4.75],...
    'ylim',[0 ceil(max(max(IntensityMeans_Norm_ROI1)))],... 
    'Box', 'off',...
    'Position', [0.08 0.5 0.25 0.25]);  % box on and off
strYLabel = '405/488';
ylabel( strYLabel, ...
    'FontName', strFontName, ...
    'FontUnit', strFontUnit, ...
    'FontSize', iFontSize, ...
    'FontWeight', strFontWeight);
title('DPI3')


% #####################
subplot_er(2,3,2)
index5dpi = [];
for jj = 1:size(Files_filtered,1) 
   if contains(Files_405_filtered(jj).name,'5dpi')
     index5dpi = [index5dpi,jj];             
   end  
end
IntensityMeans_Norm_5dpi_ROI1 = IntensityMeans_Norm_ROI1(index5dpi,:);
MeanIntensityMeans_Norm_5dpi_ROI1 = mean(IntensityMeans_Norm_5dpi_ROI1,'omitnan');
SEMIntensityMeans_Norm_5dpi_ROI1 = std(IntensityMeans_Norm_5dpi_ROI1,'omitnan')/sqrt(size(index5dpi,2));
CI95IntensityMeans_Norm_5dpi_ROI1 = [MeanIntensityMeans_Norm_5dpi_ROI1+1.96*SEMIntensityMeans_Norm_5dpi_ROI1;
    MeanIntensityMeans_Norm_5dpi_ROI1-1.96*SEMIntensityMeans_Norm_5dpi_ROI1];
% first part of figure, the individual points
% first the boxplot
bplot(IntensityMeans_Norm_5dpi_ROI1(:,1),1, 'linewidth', 4)
hold on
bplot(IntensityMeans_Norm_5dpi_ROI1(:,2),2, 'linewidth', 4)
bplot(IntensityMeans_Norm_5dpi_ROI1(:,3),3, 'linewidth', 4)
bplot(IntensityMeans_Norm_5dpi_ROI1(:,4),4, 'linewidth', 4)
% then the individual points
x_axis = repmat([1:4]',size(IntensityMeans_Norm_5dpi_ROI1,1),1);
IntensityMeans_Norm_T_5dpi_ROI1 = IntensityMeans_Norm_5dpi_ROI1(:,1:4)';
y_axis = IntensityMeans_Norm_T_5dpi_ROI1(:);
scatter(x_axis, y_axis, ...
    'MarkerEdgeColor',"#999999",...
              'MarkerFaceColor',"#999999",...
              'LineWidth',4)
line([0 4.75], [1 1],'LineWidth',3,'color',"#999999", 'LineStyle', '--')
% % finally the lines between the points
% plot(x_axis, y_axis, ':',...
%              'Color', "#999999",...
%               'LineWidth',2)
set(gca, ...
    'XTick', 1:1:4, ... ticks of x axis
    'YTick', 0:1:ceil(max(max(IntensityMeans_Norm_ROI1)))-1, ... ticks of y axis
    'XTickLabel', [],...
    'XGrid', 'off', ... [on | {off}]
    'YGrid', 'off', ... [on | {off}]
    'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorGrid', 'off' , ... [on | {off}]
    'YMinorGrid', 'off', ... [on | {off}]
    'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorTick', 'off' , ... [on | {off}]
    'YMinorTick', 'off', ... [on | {off}]
    'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
    'TickLength', [.02 .02], ... length of the ticks
    'XColor', [.1 .1 .1], ... color of x axis
    'YColor', [.1 .1 .1], ... color of y axis
    'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
    'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
    'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'FontName', strFontName, ... kind of fonts of labels
    'FontSize', iFontSize, ... size of fonts of labels
    'FontUnits', strFontUnit, ... units of the size of fonts
    'FontWeight', strFontWeight, ... weight of fonts of labels
    'FontAngle', strFontAngle, ... inclination of fonts of labels
    'LineWidth', fLineWidth, ... % width of the line of the axes
    'xlim',[0.25 4.75],...
    'ylim',[0 ceil(max(max(IntensityMeans_Norm_ROI1)))],... 
    'Box', 'off',...
    'Position', [0.38 0.5 0.25 0.25]);  % box on and off
title('DPI5')


% #############################
subplot_er(2,3,3)
index7dpi = [];
for jj = 1:size(Files_filtered,1) 
   if contains(Files_405_filtered(jj).name,'7dpi')
     index7dpi = [index7dpi,jj];             
   end  
end
IntensityMeans_Norm_7dpi_ROI1 = IntensityMeans_Norm_ROI1(index7dpi,:);
MeanIntensityMeans_Norm_7dpi_ROI1 = mean(IntensityMeans_Norm_7dpi_ROI1,'omitnan');
SEMIntensityMeans_Norm_7dpi_ROI1 = std(IntensityMeans_Norm_7dpi_ROI1,'omitnan')/sqrt(size(index7dpi,2));
CI95IntensityMeans_Norm_7dpi_ROI1 = [MeanIntensityMeans_Norm_7dpi_ROI1+1.96*SEMIntensityMeans_Norm_7dpi_ROI1;
    MeanIntensityMeans_Norm_7dpi_ROI1-1.96*SEMIntensityMeans_Norm_7dpi_ROI1];
% first part of figure, the individual points
% first the boxplot
bplot(IntensityMeans_Norm_7dpi_ROI1(:,1),1, 'linewidth', 4)
hold on
bplot(IntensityMeans_Norm_7dpi_ROI1(:,2),2, 'linewidth', 4)
bplot(IntensityMeans_Norm_7dpi_ROI1(:,3),3, 'linewidth', 4)
bplot(IntensityMeans_Norm_7dpi_ROI1(:,4),4, 'linewidth', 4)
% then the individual points
x_axis = repmat([1:4]',size(IntensityMeans_Norm_7dpi_ROI1,1),1);
IntensityMeans_Norm_T_7dpi_ROI1 = IntensityMeans_Norm_7dpi_ROI1(:,1:4)';
y_axis = IntensityMeans_Norm_T_7dpi_ROI1(:);
scatter(x_axis, y_axis, ...
    'MarkerEdgeColor',"#999999",...
              'MarkerFaceColor',"#999999",...
              'LineWidth',4)
line([0 4.75], [1 1],'LineWidth',3,'color',"#999999", 'LineStyle', '--')
% % finally the lines between the points
% plot(x_axis, y_axis, ':',...
%              'Color', "#999999",...
%               'LineWidth',2)
set(gca, ...
    'XTick', 1:1:4, ... ticks of x axis
    'YTick', 0:1:ceil(max(max(IntensityMeans_Norm_ROI1)))-1, ... ticks of y axis
    'XTickLabel', [],...
    'XGrid', 'off', ... [on | {off}]
    'YGrid', 'off', ... [on | {off}]
    'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorGrid', 'off' , ... [on | {off}]
    'YMinorGrid', 'off', ... [on | {off}]
    'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorTick', 'off' , ... [on | {off}]
    'YMinorTick', 'off', ... [on | {off}]
    'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
    'TickLength', [.02 .02], ... length of the ticks
    'XColor', [.1 .1 .1], ... color of x axis
    'YColor', [.1 .1 .1], ... color of y axis
    'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
    'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
    'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'FontName', strFontName, ... kind of fonts of labels
    'FontSize', iFontSize, ... size of fonts of labels
    'FontUnits', strFontUnit, ... units of the size of fonts
    'FontWeight', strFontWeight, ... weight of fonts of labels
    'FontAngle', strFontAngle, ... inclination of fonts of labels
    'LineWidth', fLineWidth, ... % width of the line of the axes
    'xlim',[0.25 4.75],...
    'ylim',[0 ceil(max(max(IntensityMeans_Norm_ROI1)))],... 
    'Box', 'off',...
    'Position', [0.68 0.5 0.25 0.25]);  % box on and off
title('DPI7')


% ############################
subplot_er(2,3,4)
index3dpi = [];
for jj = 1:size(Files_filtered,1) 
   if contains(Files_405_filtered(jj).name,'3dpi')
     index3dpi = [index3dpi,jj];             
   end  
end
IntensityMeans_Norm_3dpi_ROI2 = IntensityMeans_Norm_ROI2(index3dpi,:);
MeanIntensityMeans_Norm_3dpi_ROI2 = mean(IntensityMeans_Norm_3dpi_ROI2,'omitnan');
SEMIntensityMeans_Norm_3dpi_ROI2 = std(IntensityMeans_Norm_3dpi_ROI2,'omitnan')/sqrt(size(index3dpi,2));
CI95IntensityMeans_Norm_3dpi_ROI2 = [MeanIntensityMeans_Norm_3dpi_ROI2+1.96*SEMIntensityMeans_Norm_3dpi_ROI2;
    MeanIntensityMeans_Norm_3dpi_ROI2-1.96*SEMIntensityMeans_Norm_3dpi_ROI2];
% first part of figure, the individual points
% first the boxplot
bplot(IntensityMeans_Norm_3dpi_ROI2(:,1),1, 'linewidth', 4)
hold on
bplot(IntensityMeans_Norm_3dpi_ROI2(:,2),2, 'linewidth', 4)
bplot(IntensityMeans_Norm_3dpi_ROI2(:,3),3, 'linewidth', 4)
bplot(IntensityMeans_Norm_3dpi_ROI2(:,4),4, 'linewidth', 4)
% then the individual points
x_axis = repmat([1:4]',size(IntensityMeans_Norm_3dpi_ROI2,1),1);
IntensityMeans_Norm_T_3dpi_ROI2 = IntensityMeans_Norm_3dpi_ROI2(:,1:4)';
y_axis = IntensityMeans_Norm_T_3dpi_ROI2(:);
scatter(x_axis, y_axis, ...
    'MarkerEdgeColor',"#999999",...
              'MarkerFaceColor',"#999999",...
              'LineWidth',4)
line([0 4.75], [1 1],'LineWidth',3,'color',"#999999", 'LineStyle', '--')
% % finally the lines between the points
% plot(x_axis, y_axis, ':',...
%              'Color', "#999999",...
%               'LineWidth',2)
set(gca, ...
    'XTick', 1:1:4, ... ticks of x axis
    'YTick', 0:1:ceil(max(max(IntensityMeans_Norm_ROI1)))-1, ... ticks of y axis
    'XTickLabel', {'Bin1','Bin2','Bin3','Bin4'},...
    'XGrid', 'off', ... [on | {off}]
    'YGrid', 'off', ... [on | {off}]
    'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorGrid', 'off' , ... [on | {off}]
    'YMinorGrid', 'off', ... [on | {off}]
    'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorTick', 'off' , ... [on | {off}]
    'YMinorTick', 'off', ... [on | {off}]
    'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
    'TickLength', [.02 .02], ... length of the ticks
    'XColor', [.1 .1 .1], ... color of x axis
    'YColor', [.1 .1 .1], ... color of y axis
    'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
    'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
    'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'FontName', strFontName, ... kind of fonts of labels
    'FontSize', iFontSize, ... size of fonts of labels
    'FontUnits', strFontUnit, ... units of the size of fonts
    'FontWeight', strFontWeight, ... weight of fonts of labels
    'FontAngle', strFontAngle, ... inclination of fonts of labels
    'LineWidth', fLineWidth, ... % width of the line of the axes
    'xlim',[0.25 4.75],...
    'ylim',[0 ceil(max(max(IntensityMeans_Norm_ROI1)))],... 
    'Box', 'off',...
    'Position', [0.08 0.22 0.25 0.25]);  % box on and off
strYLabel = '405/488';
ylabel( strYLabel, ...
    'FontName', strFontName, ...
    'FontUnit', strFontUnit, ...
    'FontSize', iFontSize, ...
    'FontWeight', strFontWeight);


% #####################
subplot_er(2,3,5)
index5dpi = [];
for jj = 1:size(Files_filtered,1) 
   if contains(Files_405_filtered(jj).name,'5dpi')
     index5dpi = [index5dpi,jj];             
   end  
end
IntensityMeans_Norm_5dpi_ROI2 = IntensityMeans_Norm_ROI2(index5dpi,:);
MeanIntensityMeans_Norm_5dpi_ROI2 = mean(IntensityMeans_Norm_5dpi_ROI2,'omitnan');
SEMIntensityMeans_Norm_5dpi_ROI2 = std(IntensityMeans_Norm_5dpi_ROI2,'omitnan')/sqrt(size(index5dpi,2));
CI95IntensityMeans_Norm_5dpi_ROI2 = [MeanIntensityMeans_Norm_5dpi_ROI2+1.96*SEMIntensityMeans_Norm_5dpi_ROI2;
    MeanIntensityMeans_Norm_5dpi_ROI2-1.96*SEMIntensityMeans_Norm_5dpi_ROI2];
% first part of figure, the individual points
% first the boxplot
bplot(IntensityMeans_Norm_5dpi_ROI2(:,1),1, 'linewidth', 4)
hold on
bplot(IntensityMeans_Norm_5dpi_ROI2(:,2),2, 'linewidth', 4)
bplot(IntensityMeans_Norm_5dpi_ROI2(:,3),3, 'linewidth', 4)
bplot(IntensityMeans_Norm_5dpi_ROI2(:,4),4, 'linewidth', 4)
% then the individual points
x_axis = repmat([1:4]',size(IntensityMeans_Norm_5dpi_ROI2,1),1);
IntensityMeans_Norm_T_5dpi_ROI2 = IntensityMeans_Norm_5dpi_ROI2(:,1:4)';
y_axis = IntensityMeans_Norm_T_5dpi_ROI2(:);
scatter(x_axis, y_axis, ...
    'MarkerEdgeColor',"#999999",...
              'MarkerFaceColor',"#999999",...
              'LineWidth',4)
line([0 4.75], [1 1],'LineWidth',3,'color',"#999999", 'LineStyle', '--')
% % finally the lines between the points
% plot(x_axis, y_axis, ':',...
%              'Color', "#999999",...
%               'LineWidth',2)
set(gca, ...
    'XTick', 1:1:4, ... ticks of x axis
    'YTick', 0:1:ceil(max(max(IntensityMeans_Norm_ROI1)))-1, ... ticks of y axis
    'XTickLabel', {'Bin1','Bin2','Bin3','Bin4',''},...
    'XGrid', 'off', ... [on | {off}]
    'YGrid', 'off', ... [on | {off}]
    'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorGrid', 'off' , ... [on | {off}]
    'YMinorGrid', 'off', ... [on | {off}]
    'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorTick', 'off' , ... [on | {off}]
    'YMinorTick', 'off', ... [on | {off}]
    'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
    'TickLength', [.02 .02], ... length of the ticks
    'XColor', [.1 .1 .1], ... color of x axis
    'YColor', [.1 .1 .1], ... color of y axis
    'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
    'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
    'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'FontName', strFontName, ... kind of fonts of labels
    'FontSize', iFontSize, ... size of fonts of labels
    'FontUnits', strFontUnit, ... units of the size of fonts
    'FontWeight', strFontWeight, ... weight of fonts of labels
    'FontAngle', strFontAngle, ... inclination of fonts of labels
    'LineWidth', fLineWidth, ... % width of the line of the axes
    'xlim',[0.25 4.75],...
    'ylim',[0 ceil(max(max(IntensityMeans_Norm_ROI1)))],... 
    'Box', 'off',...
    'Position', [0.38 0.22 0.25 0.25]);  % box on and off



% ########################
subplot_er(2,3,6)
index7dpi = [];
for jj = 1:size(Files_filtered,1) 
   if contains(Files_405_filtered(jj).name,'7dpi')
     index7dpi = [index7dpi,jj];             
   end  
end
IntensityMeans_Norm_7dpi_ROI2 = IntensityMeans_Norm_ROI2(index7dpi,:);
MeanIntensityMeans_Norm_7dpi_ROI2 = mean(IntensityMeans_Norm_7dpi_ROI2,'omitnan');
SEMIntensityMeans_Norm_7dpi_ROI2 = std(IntensityMeans_Norm_7dpi_ROI2,'omitnan')/sqrt(size(index7dpi,2));
CI95IntensityMeans_Norm_7dpi_ROI2 = [MeanIntensityMeans_Norm_7dpi_ROI2+1.96*SEMIntensityMeans_Norm_7dpi_ROI2;
    MeanIntensityMeans_Norm_7dpi_ROI2-1.96*SEMIntensityMeans_Norm_7dpi_ROI2];
% first part of figure, the individual points
% first the boxplot
bplot(IntensityMeans_Norm_7dpi_ROI2(:,1),1, 'linewidth', 4)
hold on
bplot(IntensityMeans_Norm_7dpi_ROI2(:,2),2, 'linewidth', 4)
bplot(IntensityMeans_Norm_7dpi_ROI2(:,3),3, 'linewidth', 4)
bplot(IntensityMeans_Norm_7dpi_ROI2(:,4),4, 'linewidth', 4)
% then the individual points
x_axis = repmat([1:4]',size(IntensityMeans_Norm_7dpi_ROI2,1),1);
IntensityMeans_Norm_T_7dpi_ROI2 = IntensityMeans_Norm_7dpi_ROI2(:,1:4)';
y_axis = IntensityMeans_Norm_T_7dpi_ROI2(:);
scatter(x_axis, y_axis, ...
    'MarkerEdgeColor',"#999999",...
              'MarkerFaceColor',"#999999",...
              'LineWidth',4)
line([0 4.75], [1 1],'LineWidth',3,'color',"#999999", 'LineStyle', '--')
% % finally the lines between the points
% plot(x_axis, y_axis, ':',...
%              'Color', "#999999",...
%               'LineWidth',2)
set(gca, ...
    'XTick', 1:1:4, ... ticks of x axis
    'YTick', 0:1:ceil(max(max(IntensityMeans_Norm_ROI1)))-1, ... ticks of y axis
    'XTickLabel', {'Bin1','Bin2','Bin3','Bin4'},...
    'XGrid', 'off', ... [on | {off}]
    'YGrid', 'off', ... [on | {off}]
    'GridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorGrid', 'off' , ... [on | {off}]
    'YMinorGrid', 'off', ... [on | {off}]
    'MinorGridLineStyle', ':', ... [- | -- | {:} | -. | none]
    'XMinorTick', 'off' , ... [on | {off}]
    'YMinorTick', 'off', ... [on | {off}]
    'TickDir', 'out', ... [{in} | out] inside or outside (for 2D)
    'TickLength', [.02 .02], ... length of the ticks
    'XColor', [.1 .1 .1], ... color of x axis
    'YColor', [.1 .1 .1], ... color of y axis
    'XAxisLocation', 'bottom', ... where labels have to be printed [top | {bottom}]
    'YAxisLocation', 'left', ... where labels have to be printed [left | {right}]
    'XDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'YDir', 'normal', ... axis increasement direction [{normal} | reverse]
    'FontName', strFontName, ... kind of fonts of labels
    'FontSize', iFontSize, ... size of fonts of labels
    'FontUnits', strFontUnit, ... units of the size of fonts
    'FontWeight', strFontWeight, ... weight of fonts of labels
    'FontAngle', strFontAngle, ... inclination of fonts of labels
    'LineWidth', fLineWidth, ... % width of the line of the axes
    'xlim',[0.25 4.75],...
    'ylim',[0 ceil(max(max(IntensityMeans_Norm_ROI1)))],... 
    'Box', 'off',...
    'Position', [0.68 0.22 0.25 0.25]);  % box on and off




iResolution = 300;
% we select to crop or not the figure
set(gcf, 'PaperPositionMode', 'auto');
% saving full screen size
F    = getframe(FigH);
imwrite(F.cdata, strcat(pathpath, '_spatialTrend', '.jpeg'), 'jpeg')