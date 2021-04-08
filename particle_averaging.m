clear all; close all; clc


%% User specified parameters

% path to cropped images of indivual sites as outputed by imageJ script
dataPath = '.\data';

% wildcard identifiers to identify files of interest. If you have more than
% one experimental condition you need more than one entry
conditionStrings ={'*BRCA1*'};
%conditionStrings ={'**/*Mid S*/*', '**/*Late S*/*'};

% condition labels for plotting
plotStrings ={'BRCA1'};
%plotStrings ={'Mid S', 'Late S'};

% dimensions of voxels (in microns)
voxelXY = 0.159;
voxelZ = 0.2;

% bin edges for radial plots
binEdges = 0:0.2:2.4;

% classes to process (as defined by imageJ script
classes = [5];


%% Main script
        
% retrieve working directory containing the script
workDir = pwd;

% number of conditions, classes and plotting bins as specified by user
numClasses = length(classes);
numConditions = length(conditionStrings);
numBins = length(binEdges) - 1;

% create blank figure windows for plots
radFig = figure;
c1Plot = figure;
c2Plot = figure;

% loop through conditions and classes
for cn = 1 : numConditions

    for cl = 1 : numClasses

        % go to directory containing data
        cd(dataPath);
        % find all site crops correspoding to the condition and class
        pathCard = strcat(conditionStrings{cn},'class', num2str(classes(cl)), '*.tif'); 
        files = dir(pathCard);
        % go back to working directory
        cd(workDir);
        
        % number of sites to process for this condition and class
        numFiles = size(files, 1);
        
        % stack to hold cropped sites
        particleStack = [];
        % count of sites loaded
        count = 0;
        for i = 1 : numFiles
            % load site image
            filePath = strcat(files(i).folder, filesep, files(i).name);  
            fileInfo = imfinfo(filePath);
            numImages = numel(fileInfo);
            % if site image width = height, ie not only image edge
            if (fileInfo(1).Width == fileInfo(1).Height)
                count = count + 1;
                channel = 1;
                slice = 1;
                for k = 1 : numImages
                    particleStack(:, :, channel, slice, count) = imread(filePath, k);
                    if channel == 1
                        channel = 2;
                    else
                        slice = slice + 1;
                        channel = 1;
                    end
                end 
            end

        end
        % compute average of all sites
        particleAverage = mean(particleStack, 5);
        
        % retrieve dimensions
        [height, width, numChannels, numSlices] = size(particleAverage);
    
        % create coordinate grids for X and Y centred on site (scaled to microns)
       [X, Y] = meshgrid((1 : width) - ceil(width / 2),(1 : height) - ceil(height / 2));
        X = X * voxelXY;
        Y = Y * voxelXY;
        
        % radial distance from site centre (scaled to microns) 
        R = sqrt(X.^2 + Y.^2);
        
        % loop through both channels
        for c = 1 : 2
            
            % to contain radial profile and bin centres
            radialProfile = [];
            binCentres = [];

            % loop through bins
            for b = 1 : numBins
               % mask for current bin
               mask = (binEdges(b) <= R & R < binEdges(b + 1));
               % central slice of averaged site
               band = particleAverage(:, :, c, ceil(numSlices / 2));
               % vectorise both
               mask = mask(:);
               band = band(:);
               % delete entries of band outside of mask
               band(mask == false) = [];
               % retrieve mean intensity of band 
               radialProfile(b) = mean(band);
               % centre of bin (for plotting)
               binCentres(b) = (binEdges(b + 1) - binEdges(b)) /2 + binEdges(b);
            end

            % select radial plots
            figure(radFig)
            % create/select subplot for this condition, class and channel
            subplot(numConditions, 2, (cn - 1) * 2 + c);
            hold on
            % plot radial profile
            plot(binCentres, radialProfile, 'LineWidth', 2);
            % if last class add a legend and title
            if (cl == numClasses) 
                legend(string(classes))
                xlabel('radius (microns)')
                ylabel('mean intensity')
                title(strcat(plotStrings{cn}, ' channel ', num2str(c)))
                % axis plot limits (different for each channel)
                if (c == 1)
                    ylim([0 8])
                end
                if (c == 2)
                    ylim([0 40])
                end
            end
            
            % select figure for particle average corresponding to current
            % channel
            if (c == 1)
                figure(c1Plot)
            else
                figure(c2Plot)
            end
            % create/select subplot for this condition, class and channel
            subplot(numConditions, numClasses, (cn - 1) * (numClasses) + cl );
            
            % shows particle average at central slic
            imshow(particleAverage(:, :, c, ceil(numSlices / 2)), [])
            % add sugplot title
            plotTitle = strcat('class', num2str(classes(cl)), ', channel ', num2str(c), ', ', plotStrings{cn});
            title(plotTitle)
            
            % calculate min and max value of particle average for channel
            channelAverage = particleAverage(:,:,c,:);
            channelMax = max(channelAverage(:));
            channelMin = min(channelAverage(:));
            % loop though all slices
            for k=1:numSlices     
                % normalise slice
                sliceNorm = (particleAverage(:,:,c,k) - channelMin) / (channelMax - channelMin);
                % save particle average to tif stack
                if k == 1
                    imwrite(sliceNorm, strcat(plotTitle, '.tif'));    
                else
                    imwrite(sliceNorm, strcat(plotTitle, '.tif') ,'WriteMode','append');    
                end
                
            end
        end
    end

 


end