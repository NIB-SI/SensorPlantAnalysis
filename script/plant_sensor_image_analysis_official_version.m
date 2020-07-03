%% 
% plant_sensor_image_analysis_official_version.m
% Date: Jan 29, 2020
% Author: Anze Zupanic
% Responsible for test data: Tjasa Lukan
%                 these data are from 2019, so the comments are specific to
%                 these data only. Ask Tjasa for more info on how they were
%                 obtained. 

%% USER IMPORTANT INFORMATION: 5 stars ***** means that the code needs adjustment (correct path, filenames, etc) !!!!!!!

%% IMPORT AND DETECTION OF CHLOROPLASTS
    % extracts the chloroplast and creates a mask
    % Output Image: image diary of the analysis steps taken
    % Output Variables:
        % masks - mask image itself
        % chloroplastArea - the are of chloroplasts in the image 
        % pixels - pixels belonging to each chloroplast (I_connected)
        % fluorescence - fluorescent intensity of each identified chloroplast (grayscaleValues)

% ***** MAKE SURE THAT THE PATH IS CORRECT ******
pathpath = 'O:\DEJAVNOSTI\OMIKE\Methods_DryLab\Task_ImageAnalysis\Image_analysis_plant_sensors\testData\exp1_chlorophyll_originals\';
% reads all tiff files (= chloroplast images)
Files = dir(strcat(pathpath,'*.tif'));
numFiles = numel(Files); % returns number of files found

% loop over all files
for jj = 1:numFiles

    jj
    % if at any point you want to plot the image, use imshow(image_name)
   
    % read in each image and convert to grayscale
    image_original = imread(strcat(Files(jj).folder,'\',Files(jj).name));
    image_gray = rgb2gray(image_original);
    
    % first a very simplistic way of getting rid of noise (background signal
    % from other cells,etc)
    % figures are reasonably empty, but really noisy, so it-s not possible
    % to learn about signal vs background from the histogram alone, but the
    % peak is usually in background, so this code get's rid of much of the
    % noise, without getting rid of any of the signal
    [counts,x] = imhist(image_gray);
    [i,whereMax] = max(counts);
    % ***** try other values here
    image_gray(image_gray<whereMax+5)=0;
    
    
    % make a binary image, based on a local adaptive method
    %   sensitivity = 0.6 -
    %   neighbourhood size = 15 (will  separate between the
    %   chloroplasts)
    % the idea here is to include as much of the signal as possible, and
    % remove almost all of the noise
    % ***** (try playing with settings a bit)
    T = adaptthresh(image_gray,0.6,'NeighborhoodSize',11, 'Statistic', 'Gaussian');
    image_binary = imbinarize(image_gray, T);
    % old code, which only looked at sensitivity *just for comparison)
    % image_binary = imbinarize(image_gray, 'adaptive');
    
    % erosion, dilation, median filtering and areaopening operations provide
    % good results, regardless of image quality
    se = strel('square', 6);
    image_open = imopen(image_binary,se);
    image_filtered = medfilt2(image_open);
    % here the minimum size of the chloroplast is defined, anything smaller
    % gets deleted
    % *****define the minimum size according to species/imaging protocol
    image_sized = bwareaopen(image_filtered, 120);
    % so, now we have a mask that has almost individual chloroplasts, and 
    % and is losing almost no deep chloroplasts. for something better,
    % machine learning would probably have to be used
    
    % OUTPUT - masks
    masks{jj} = image_sized;
    
    % calculate total chloroplast are in the image (area of mask)
    chloroplastArea{jj} = bwarea(image_sized);

    % NORMALIZATION
    % calculate mean fluorescence of chloroplasts
    %   - change image to double * 255 gives the same values as in image_gray
    %   - then multiply with mask, now you have 0s where no mask and gray 
    %       values where mask
    %   - mean(mean*image_temp)) now gives the average value of the whole
    %       image (and not only mask), so it needs to be normalized 
    %       ***** 1048576 is total number of pixels in image (change this if zour images are of different size)
    %       chloroplastArea is total number of pixels in mask
    image_temp = im2double(image_gray).*255;
    image_temp = image_temp.*image_sized;
    meanChloroplastFluorescence{jj} = 1048576*mean(mean(image_temp))/chloroplastArea{jj};
    
    % now we have the mask, it-s time to get the individual chloroplasts out,
    % plus clusters of chloroplasts where it-s difficult to divide them into
    % individual ones.
    % chloroplast_pixels.PixelldxList - indexes where the pixels belong to
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
    % contour plot, i don-t like it, but it stays here for future use
%     subplot(3,3,7)
%     imcontour(image_sized)
%     title('contour plot')
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
    to show the figure
    figure(FigH)
    
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

%saving the important results
savefile = strcat(pathpath,'analyzedChloroplasts.mat');
% here we define all OUTPUTs to save
save(savefile, 'masks','chloroplast_pixels','chloroplastStats', 'figureStats', 'Files');

% loading the important results
% ***** define your path
path = 'O:\DEJAVNOSTI\OMIKE\Methods_DryLab\Task_ImageAnalysis\Image_analysis_plant_sensors\testData\exp1_chlorophyll_originals\';
load(strcat(path,'analyzedChloroplasts.mat'))


%% 405/488 RATIO per FIGURE and per CHLOROPLAST
    % masks from first part are used on the 405 and 488 datasets
    % the ratio between the two fluorescences is calculated
        
% get all the images of 405 and 488
% it needs to be the same number of images as for the chlorophyll
% ***** MAKE SURE THAT THE PATH IS CORRECT *****   
Files_405 = dir('O:\DEJAVNOSTI\OMIKE\Methods_DryLab\Task_ImageAnalysis\Image_analysis_plant_sensors\testData\exp1_fluo405_originals\*.tif');
Files_488 = dir('O:\DEJAVNOSTI\OMIKE\Methods_DryLab\Task_ImageAnalysis\Image_analysis_plant_sensors\testData\exp1_fluo488_originals\*.tif');
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
    
    % calculate the ratio between the images, without any changes 405/488
    % division comes first, then the application of mask
    image_rat = image_405_double./image_488_double;
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




%% 405/488 RATIO per FIGURE and per CHLOROPLAST, normalized against the background
    % masks from first part are used on the 405 and 488 datasets
    % the ratio between the two fluorescences is calculated
        
% get all the images of 405 and 488
% it needs to be the same number of images as for the chlorophyll
% ***** MAKE SURE THAT THE PATH IS CORRECT *****   
Files_405 = dir('O:\DEJAVNOSTI\OMIKE\Methods_DryLab\Task_ImageAnalysis\Image_analysis_plant_sensors\testData\exp1_fluo405_originals\*.tif');
Files_488 = dir('O:\DEJAVNOSTI\OMIKE\Methods_DryLab\Task_ImageAnalysis\Image_analysis_plant_sensors\testData\exp1_fluo488_originals\*.tif');
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


