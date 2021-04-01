clear all; close all; clc


%% User specified parameters
dataPath = 'D:\SMARCAD1 depletions U2OS empty\Class 5 foci\NTC\Normal exclusion\';


% dataPath = 'E:\Visualising BRCA1 and 53BP1 foci 3 colour\';


conditionStrings ={'**/*Mid S*/*', '**/*Late S*/*'};


%conditionStrings ={'*BRCA1*'}
%plotStrings ={'test'};
plotStrings ={'Mid S', 'Late S'};

voxelXY = 0.159;
voxelZ = 0.2;


binEdges = 0:0.2:2.4;

classes = [5];
numClasses = length(classes);


%% Main script
numConditions = length(conditionStrings);

radFig = figure;

c1Plot = figure;



c2Plot = figure;

for cn = 1 : numConditions
    
   

    for cl = 1 : numClasses
        workDir = pwd;
        cd(dataPath);
        pathCard = strcat(conditionStrings{cn},'class', num2str(classes(cl)), '*.tif'); 
        files = dir(pathCard);
        cd(workDir);

        numFiles = size(files, 1);

        particleStack = [];

        count = 0;
        for i = 1 : numFiles
            filePath = strcat(files(i).folder, filesep, files(i).name);  
            fileInfo = imfinfo(filePath);
            numImages = numel(fileInfo);
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
        
        particleAverage = mean(particleStack, 5);

        [height, width, numChannels, numSlices] = size(particleAverage);
        
        
        
        
       % [X, Y, Z] = meshgrid((1 : width) - ceil(width / 2),(1 : height) - ceil(height / 2), (1 : numSlices) - ceil(numSlices / 2));
        [X, Y] = meshgrid((1 : width) - ceil(width / 2),(1 : height) - ceil(height / 2));
        X = X * voxelXY;
        Y = Y * voxelXY;
       %  Z = Z * voxelZ;
        % R = sqrt(X.^2 + Y.^2 + Z.^2);
        R = sqrt(X.^2 + Y.^2);
        for c = 1 : 2


            radialProfile = [];
            binCentres = [];

            numBins = length(binEdges) - 1;
            for b = 1 : numBins
               mask = (binEdges(b) <= R & R < binEdges(b + 1));
               mask = mask(:);
               band = particleAverage(:, :, c, ceil(numSlices / 2)); 
               band = band(:);
               band(mask == false) = [];



               radialProfile(b) = mean(band);
               binCentres(b) = (binEdges(b + 1) - binEdges(b)) /2 + binEdges(b);
               %bandAreas(b) = pi * (binEdges(b + 1)^2 - binEdges(b)^2 );
            end

            %radialProfileNorm = radialProfile ./ bandAreas;
            
            figure(radFig)
          
            subplot(numConditions, 2, (cn - 1) * 2 + c);
            hold on
            plot(binCentres, radialProfile, 'LineWidth', 2);
            if (cl == numClasses) 
                legend(string(classes))
                xlabel('radius (microns)')
                ylabel('mean intensity')
                title(strcat(plotStrings{cn}, ' channel ', num2str(c)))
                if (c == 1)
                    ylim([0 8])
                end
                   if (c == 2)
                    ylim([0 40])
                end
            end

        %     figure; 
        %     plot(binCentres, radialProfileNorm);
        %     title(['norm class' num2str(cl)])
            if (c == 1)
                figure(c1Plot)
            else
                figure(c2Plot)
            end
            (cn - 1) * numConditions + cl
            subplot(numConditions, numClasses, (cn - 1) * (numClasses) + cl );
            
            
            % shows particle image 
            imshow(particleAverage(:, :, c, ceil(numSlices / 2)), [])
            plotTitle = strcat('class', num2str(classes(cl)), ', channel ', num2str(c), ', ', plotStrings{cn});
            title(plotTitle)
            for k=1:numSlices     
                channelAverage = particleAverage(:,:,c,:);
                channelMax = max(channelAverage(:));
                channelMin = min(channelAverage(:));
                sliceNorm = (particleAverage(:,:,c,k) - channelMin) / (channelMax - channelMin);
                if k == 1
                    imwrite(sliceNorm, strcat(plotTitle, '.tif'));    
                else
                    imwrite(sliceNorm, strcat(plotTitle, '.tif') ,'WriteMode','append');    
                end
                
            end
        end
    end

 


end