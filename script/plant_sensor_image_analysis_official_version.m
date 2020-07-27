%% 
% plant_sensor_image_analysis_official_version.m
% Date of last change: July 17, 2020
% Author: Anze Zupanic
% Responsible for test data: Tjasa Lukan
%                 these data are from 2019, so the comments are specific to
%                 these data only. Ask Tjasa for more info on how they were
%                 obtained.

%% USER IMPORTANT INFORMATION: 5 stars ***** means that the code needs adjustment (correct path, filenames, etc) !!!!!!!

%% DATA AND IMAGE PREPARATION FOR THE ANALYSIS
% to perform all the different steps of the analysis in this script, the
% following is required:
    % microscopic images, 3 different channels, chlorophyl channel, 405 nm
    % channel and 408 nm channel. IMPORTANT: there should be a separate
    % folder for each type of images, and the names of the files should be
    % the SAME, except for the ending of the files, these are good names: 
        % 17022020_Rywal_NahG_pt-roGFP_L2_P1_L2_PVY-WILGA_3dpi.lif_lesion01_ROI01_Processed001_ch00.tif
        % 17022020_Rywal_NahG_pt-roGFP_L2_P1_L2_PVY-WILGA_3dpi.lif_lesion01_ROI01_Processed001_ch01.tif
        % 17022020_Rywal_NahG_pt-roGFP_L2_P1_L2_PVY-WILGA_3dpi.lif_lesion01_ROI01_Processed001_ch02.tif
    % for the spatial analysis, each file name (only one channel can be
    % used here, it's just important that the alphabetical order matches),
    % should be accompanied by the position of the lesion
% example files, with appropriate naming and formats can be found on 
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
        
% ***** MAKE SURE THAT THE PATH IS CORRECT *****
% path to the chloroplast files
pathpath = 'C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch00_chlorophyll\';
% reads all tiff files (= chloroplast images)
Files = dir(strcat(pathpath,'*.tif'));
numFiles = numel(Files); % returns number of files found

% loop over all files
for jj = 1:numFiles
     
    % read in each image and convert to grayscale
    image_original = imread(strcat(Files(jj).folder,'\',Files(jj).name));
    image_gray = rgb2gray(image_original);
    
    % first a very simplistic way of getting rid of noise (background signal
    % from other cells,etc). Finds the peak of the histogram and deleted
    % everything to the right of the peak +some offset
    [counts,x] = imhist(image_gray);
    [i,whereMax] = max(counts);
    % ***** try other values here
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
    image_sized = bwareaopen(image_filtered, 120);
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
    %       ***** 1048576 is total number of pixels in image (change this if your images are of different size)
    image_temp = im2double(image_gray).*255;
    image_temp = image_temp.*image_sized;
    meanChloroplastFluorescence{jj} = 1048576*mean(mean(image_temp))/chloroplastArea{jj};
    
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
    %***** change this to your saving folder*****
    strFilePath = 'H:\projekti\projekti NIB\lezije\exp1_chlorophyll_originals';
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


% loading the important results
% ***** define your path
pathpath = 'C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch00_chlorophyll\';
load(strcat(pathpath,'analyzedChloroplasts.mat'))


%% PART 2A: 405/488 RATIO per FIGURE and per CHLOROPLAST
    % masks from first part are used on the 405 and 488 microscopic images
    % the ratio between the two fluorescences is calculated
        
% get all the images of 405 and 488
% it needs to be the same number of images as for the chlorophyll
% ***** MAKE SURE THAT THE PATH IS CORRECT *****   
pathpath = 'C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch00_chlorophyll\';
Files_405 = dir('C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch01_405\*.tif');
Files_488 = dir('C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch02_488\*.tif');
load(strcat(pathpath,'analyzedChloroplasts.mat'))
numFiles = numel(Files_405);
% read the images one by one and do the analysis
for jj = 1:numFiles
    
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
    % ***** 1048576 is total number of pixels in image (change this if your images are of different size)
    normalizedRatio{jj} = 1048576*mean(mean(image_ratio))/figureStats.AreaChloroplasts(jj);
    
    
%     % PLOTTING
%     % series of images to illustrate, what was done
%     % ***** uncomment this if you want plotting and saving the images
%     FigH = figure('Position', get(0, 'Screensize'));
%     set(FigH, 'Visible', 'off');
%     subplot(2,2,1)
%     imshow(image_405_gray)
%     caxis([0 30])
%     title('Scaled 405')
%     subplot(2,2,2)
%     imshow(image_488_gray)
%     caxis([0 30])
%     title('Scaled 488')
%     subplot(2,2,3)
%     image(image_rat,'CDataMapping','scaled')
%     axis equal
%     colormap(gray)
%     % ***** check the histogram and then define the color axis
%     % hist(image_rat(:),100)
%     caxis([0 4])
%     title('405/488')
%     subplot(2,2,4)
%     image(image_ratio,'CDataMapping','scaled')
%     axis equal
%     colormap(gray)
%     % ***** check the histogram and then define the color axis
%     % hist(image_ratio(:),100)
%     caxis([0 4])
%     title('405/488 in chloroplasts')
%     % to show the figure
%     figure(FigH)
    
%     % select path and save
%     %***** change this to your saving folder*****
%     strFilePath = 'H:\projekti\projekti NIB\lezije\exp1_chlorophyll_originals';
%     % we select the printing resolution
%     iResolution = 300;
%     % we select to crop or not the figure
%     set(gcf, 'PaperPositionMode', 'auto');
%     % saving full screen size
%     F    = getframe(FigH);
%     imwrite(F.cdata, strcat(pathpath, Files(jj).name, '_ratio','.jpeg'), 'jpeg')
     

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




%% PART 2B: 405/488 RATIO per FIGURE and per CHLOROPLAST, NORMALIZED AGAINST BACKGROUND
    % the idea here is that the signal in 488 is so low, that even small
    % changes in background between images can bias the analysis. The code
    % here is the same as above in PART 2A, but here an average background
    % in 488nm is found, and then removed from the image. However, this
    % average background can depend on the number of chloroplasts in the
    % image (there is a halo around each chloroplast), so in the future
    % manual determination of the background would be a better solution.
    % However, if this simple the removal of background does not actually 
    % have a significant impact, then manual determination is a waste of
    % time. 
        
% get all the images of 405 and 488
% it needs to be the same number of images as for the chlorophyll
% ***** MAKE SURE THAT THE PATH IS CORRECT *****   
pathpath = 'C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch00_chlorophyll\';
Files_405 = dir('C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch01_405\*.tif');
Files_488 = dir('C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch02_488\*.tif');
load(strcat(pathpath,'analyzedChloroplasts.mat'))
numFiles = numel(Files_405);
% read the images one by one and do the analysis
for jj = 1:numFiles
    jj
    
    % PER FIGURE
    % read in the image and convert to grayscale
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
    
    
    % HERE NORMALIZATION AGAINST BACKGROUND FOR 488
    % calculate the background 488 for each image
    image_488_background = image_488_double.*imcomplement(masks{jj});
    average_background = 1048576*mean(mean(image_488_background))/(1048576 - figureStats.AreaChloroplasts(jj));
    image_488_removedBackground = image_488_double - average_background;
    % to prevent negative elements and values smaller than 1
    image_488_removedBackground(image_488_removedBackground < 1) = 1;
    
    % calculate the ratio between the images, without any changes 405/488
    % division comes first, then the application of mask
    image_rat = image_405_double./image_488_removedBackground;
    image_ratio = image_rat.*masks{jj};
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
    % ***** 1048576 is total number of pixels in image (change this if zour images are of different size)
    normalizedRatio_background{jj} = 1048576*mean(mean(image_ratio))/figureStats.AreaChloroplasts(jj);
    
    
%     % PLOTTING
%     % series of images to illustrate, what was done
%     % ***** uncomment this if you want plotting and saving the images
%     FigH = figure('Position', get(0, 'Screensize'));
%     set(FigH, 'Visible', 'off');
%     subplot(2,2,1)
%     imshow(image_405_gray)
%     caxis([0 30])
%     title('Scaled 405')
%     subplot(2,2,2)
%     imshow(image_488_gray)
%     caxis([0 30])
%     title('Scaled 488')
%     subplot(2,2,3)
%     image(image_rat,'CDataMapping','scaled')
%     axis equal
%     colormap(gray)
%     % ***** check the histogram and then define the color axis
%     % hist(image_rat(:),100)
%     caxis([0 4])
%     title('405/488')
%     subplot(2,2,4)
%     image(image_ratio,'CDataMapping','scaled')
%     axis equal
%     colormap(gray)
%     % ***** check the histogram and then define the color axis
%     % hist(image_ratio(:),100)
%     caxis([0 4])
%     title('405/488 in chloroplasts')
%     % to show the figure
%     figure(FigH)
    
%     % select path and save
%     %***** change this to your saving folder*****
%     strFilePath = 'H:\projekti\projekti NIB\lezije\exp1_chlorophyll_originals';
%     % we select the printing resolution
%     iResolution = 300;
%     % we select to crop or not the figure
%     set(gcf, 'PaperPositionMode', 'auto');
%     % saving full screen size
%     F    = getframe(FigH);
%     imwrite(F.cdata, strcat(pathpath, Files(jj).name, '_ratio','.jpeg'), 'jpeg')
     

    % PER CHLOROPLAST
    % OUTPUT
    % her the ratio per all chloroplasts per each figure are saved
    grayscaleChloroplasts = regionprops(masks{jj},image_ratio,{'Centroid', 'MeanIntensity'});
    chloroplast_ratios_background{jj} = grayscaleChloroplasts;

end

% PLOTTING the results PER FIGURE
% it is very important that the figures that are analyzed are in a sensible
% order, because the points on the plot are in the same order as the
% figures that were analyzed!
normRatio = cell2mat(normalizedRatio_background);
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
   temp = [chloroplast_ratios_background{jj}.MeanIntensity];
   bplot(temp,jj)
end
set(gca,'FontSize',20)
hold off

%% PART 3: SPATIAL ANALYSIS  INCLUDED FOR THE REGIONS OF INTEREST CLOSEST TO THE LESION (ROI1), NORMALIZED AGAINST background
% this part can be performed only after finding the chloroplasts
% which means that you first need to define paths to files, and then load
% the saved data from chloroplast analysis
clear all
pathpath = 'C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch00_chlorophyll\';
Files_405 = dir('C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch01_405\*.tif');
Files_488 = dir('C:\projekti\projekti NIB\lezije\spatialAnalysis\Export\L2\ch02_488\*.tif');
load(strcat(pathpath,'analyzedChloroplasts.mat'))
numFiles = numel(Files_405);

% find all images, which have ROI1 in their name and not MOCK, DTT, h202 or
% control
index = [];
for jj = 1:size(Files,1) 
   if contains(Files_405(jj).name,'ROI01') && ~contains(Files_405(jj).name,{'MOCK','DTT','H2O2','control'})
   % here for ROI2
 % if contains(Files_405(jj).name,'ROI02') && ~contains(Files_405(jj).name,{'MOCK','DTT','H2O2','control'})
     index = [index,jj];             
   end  
end
Files = Files(index,:);
Files_405 = Files_405(index,:);
Files_488 = Files_488(index,:);
chloroplastStats = chloroplastStats(1,index);
figureStats = figureStats(index,:);
chloroplast_pixels = chloroplast_pixels(index);
masks = masks(index);

% load the excel file and read in the positions of each lesion for each ROI01
% the files in the excel should be named the same as the files in the folders, but without the ending Processed_ch00 etc 
% ***** path to your excel file
positions = readtable('C:\projekti\projekti NIB\lezije\spatialAnalysis\RywalNahG_2nd_exp_lesion_position.xlsx');
% we will now add the distance to the lesion to the chloroplastStats
% variable
% ***** size of the image in this case is 512 pixels, for images of different
% size, write a different number
% for each image
% we repeat the same analysis as before, but with spatial information
for jj = 1:size(Files,1)
    % in this first part the intensities per chloroplast are calculated as
    % before
    
    % PER FIGURE
    % read in the image and convert to grayscale
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
    
    
    % HERE NORMALIZATION AGAINST BACKGROUND FOR 488
    % ***** if you don't want normalization, just comment our the next 6
    % lines
    % calculate the background 488 for each image
    image_488_background = image_488_double.*imcomplement(masks{jj});
    average_background = 1048576*mean(mean(image_488_background))/(1048576 - figureStats.AreaChloroplasts(jj));
    image_488_removedBackground = image_488_double - average_background;
    % to prevent negative elements and values smaller than 1
    image_488_removedBackground(image_488_removedBackground < 1) = 1;
    
    
    % calculate the ratio between the images, without any changes 405/488
    % division comes first, then the application of mask
    image_rat = image_405_double./image_488_removedBackground;
    image_ratio = image_rat.*masks{jj};
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
    % ***** 1048576 is total number of pixels in image (change this if zour images are of different size)
    normalizedRatio_background{jj} = 1048576*mean(mean(image_ratio))/figureStats.AreaChloroplasts(jj);
    
%     % PLOTTING
%     % series of images to illustrate, what was done
%     % ***** uncomment this if you want plotting and saving the images
%     FigH = figure('Position', get(0, 'Screensize'));
%     set(FigH, 'Visible', 'off');
%     subplot(2,2,1)
%     imshow(image_405_gray)
%     caxis([0 30])
%     title('Scaled 405')
%     subplot(2,2,2)
%     imshow(image_488_gray)
%     caxis([0 30])
%     title('Scaled 488')
%     subplot(2,2,3)
%     image(image_rat,'CDataMapping','scaled')
%     axis equal
%     colormap(gray)
%     % ***** check the histogram and then define the color axis
%     % hist(image_rat(:),100)
%     caxis([0 4])
%     title('405/488')
%     subplot(2,2,4)
%     image(image_ratio,'CDataMapping','scaled')
%     axis equal
%     colormap(gray)
%     % ***** check the histogram and then define the color axis
%     % hist(image_ratio(:),100)
%     caxis([0 4])
%     title('405/488 in chloroplasts')
%     % to show the figure
%     figure(FigH)
    
%     % select path and save
%     %***** change this to your saving folder*****
%     strFilePath = 'H:\projekti\projekti NIB\lezije\exp1_chlorophyll_originals';
%     % we select the printing resolution
%     iResolution = 300;
%     % we select to crop or not the figure
%     set(gcf, 'PaperPositionMode', 'auto');
%     % saving full screen size
%     F    = getframe(FigH);
%     imwrite(F.cdata, strcat(pathpath, Files(jj).name, '_ratio','.jpeg'), 'jpeg')
     
    % PER CHLOROPLAST
    % OUTPUT
    % her the ratio per all chloroplasts per each figure are saved
    grayscaleChloroplasts = regionprops(masks{jj},image_ratio,{'Centroid', 'MeanIntensity'});
    chloroplast_ratios_background{jj} = grayscaleChloroplasts;
    
    % in this second part, the distance between each chloroplast and the
    % lesion is calculated
    
    numChloroplasts = size(chloroplastStats{1,jj});
    % for each chloroplast
    for kk = 1:numChloroplasts(1)
        coordinatesChloroplast = chloroplastStats{1,jj}(kk).Centroid;
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

% PLOTTING the results PER FIGURE
% it is very important that the figures that are analyzed are in a sensible
% order, because the points on the plot are in the same order as the
% figures that were analyzed!
normRatio = cell2mat(normalizedRatio_background);
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
IntensityMeans = zeros(size(Files,1),5);
IntensityVariability = zeros(size(Files,1),5);
for jj = 1:size(Files,1)
    tempStruct = chloroplast_ratios_background{jj};
    % this is a very important line of code that show how to extract a
    % field from s structure (inside parentheses)
    tempVector = [tempStruct.Bin];
    for kk = 1:5 % **** because five bins
       index = tempVector == kk; 
       temp = [tempStruct.MeanIntensity];
       temp = temp(index);
       IntensityMeans(jj,kk)=mean(temp);
       IntensityVariability(jj,kk)=var(temp); 
       bplot(temp,(jj-1)*5+kk) %again 5 bins *****
    end
end
set(gca,'FontSize',20)
hold off

% mean, variability and CI for fluorescence intensities (ratios)
MeanIntensityMeans = mean(IntensityMeans,'omitnan');
SEMIntensityMeans = std(IntensityMeans,'omitnan')/sqrt(size(Files,1));
CI95IntensityMeans = [MeanIntensityMeans+1.96*SEMIntensityMeans;
    MeanIntensityMeans-1.96*SEMIntensityMeans];

% same, but normalized to BIN1
MeanIntensityMeans_norm = MeanIntensityMeans/MeanIntensityMeans(1);
SEMIntensityMeans_norm = SEMIntensityMeans/MeanIntensityMeans(1);
CI95IntensityMeans_norm = CI95IntensityMeans/MeanIntensityMeans(1);
