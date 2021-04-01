clear all; close all; clc


%% User specified parameters
dataPath = 'D:\SMARCAD1 depletions H2A Dendra\Class 5 foci structures\NTC\53BP1_BRCA1\';
addpath(genpath('receptor-trafficking-toolbox'))

% dataPath = 'E:\Visualising BRCA1 and 53BP1 foci 3 colour\';


conditionStrings ={'**/*Mid S*/*', '**/*Late S*/*'};


% conditionStrings ={'**/Mid S*/*BRCA1*', '**/Late S*/*BRCA1*'}
plotStrings ={'Mid S', 'Late S'};

numChannels = 2;
colocChannel = [1 2];

voxelXY = 0.159;
voxelZ = 0.2;


binEdges = 0:0.1:2.4;

classes = [5];
numClasses = length(classes);


%% Main script
numConditions = length(conditionStrings);

radFig = figure;

cPlots = [];
for c = 1 : numChannels
    cPlots = [cPlots figure];
end

colocStatsMean = [];
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
            
                slice = 1;
                channel = 1;
                for k = 1 : numImages
          
                        
                    particleStack(:, :, slice, 1, channel, count) = imread(filePath, k);
                    channel = channel + 1;
                    if (channel > numChannels)
                        slice = slice + 1;
                        channel = 1;
                    end
                 
                    
                end 
            end

        end
        numParticles = count;
        particleAverage = mean(particleStack, 6);

        [height, width, numSlices, numTimePoints, numChannels] = size(particleAverage);
        
        
        colocStats = [];
         for i = 1 : numParticles
             particle = particleStack(:,:,:,:,:,i);
             thresholds = otsuThresholds(particle, ones(height, width, numSlices, numTimePoints));
             [ROI, ~] = bwconvhull3d(particle(:,:,:,:,colocChannel(1)) >= thresholds(1) | particle(:,:,:,:,colocChannel(1)) >= thresholds(2)); 

             
             colocStats = [colocStats; struct2table(calcColocStats(particle, ROI, thresholds, true))];
         end
       
        colocStatsMeanTemp = [table(conditionStrings(cn), classes(cl), numParticles, 'VariableNames',{'condition', 'class', 'number'}) varfun(@nanmean, colocStats) varfun(@nanstd, colocStats)];
         colocStatsMean =[colocStatsMean; colocStatsMeanTemp];
       % [X, Y, Z] = meshgrid((1 : width) - ceil(width / 2),(1 : height) - ceil(height / 2), (1 : numSlices) - ceil(numSlices / 2));
        [X, Y] = meshgrid((1 : width) - ceil(width / 2),(1 : height) - ceil(height / 2));
        X = X * voxelXY;
        Y = Y * voxelXY;
       %  Z = Z * voxelZ;
        % R = sqrt(X.^2 + Y.^2 + Z.^2);
        R = sqrt(X.^2 + Y.^2);
        for c = 1 : numChannels


            radialProfile = [];
            binCentres = [];

            numBins = length(binEdges) - 1;
            for b = 1 : numBins
               mask = (binEdges(b) <= R & R < binEdges(b + 1));
               mask = mask(:);
               band = particleAverage(:, :, ceil(numSlices / 2), :, c); 
               band = band(:);
               band(mask == false) = [];



               radialProfile(b) = mean(band);
               binCentres(b) = (binEdges(b + 1) - binEdges(b)) /2 + binEdges(b);
               %bandAreas(b) = pi * (binEdges(b + 1)^2 - binEdges(b)^2 );
            end
            radialProfile = (radialProfile - min(radialProfile(:)))/(max(radialProfile(:))-min(radialProfile(:)))
            %radialProfileNorm = radialProfile ./ bandAreas;
            
            figure(radFig)
          
            subplot(numConditions, 1, cn);
            hold on
            plot(binCentres, radialProfile, 'LineWidth', 2);

            if (cl == numClasses) 
                
                xlabel('radius (microns)')
                ylabel('mean intensity (normalised)')
                title(plotStrings{cn})
%                 if (c == 1)
%                     ylim([0 8])
%                 end
%                    if (c == 2)
%                     ylim([0 60])
%                 end
                ylim([-0.2 1.2])
            end

        %     figure; 
        %     plot(binCentres, radialProfileNorm);
        %     title(['norm class' num2str(cl)])
           figure(cPlots(c))
            (cn - 1) * numConditions + cl
            subplot(numConditions, numClasses, (cn - 1) * (numClasses) + cl );
            
            
            % shows particle image 
            imshow(particleAverage(:, :, ceil(numSlices / 2), 1, c), [])
            plotTitle = strcat('class', num2str(classes(cl)), ', channel ', num2str(c), ', ', plotStrings{cn});
            title(plotTitle)
            for k=1:numSlices     
                channelAverage = particleAverage(:,:,:,1,c);
                channelMax = max(channelAverage(:));
                channelMin = min(channelAverage(:));
                sliceNorm = (particleAverage(:,:,k,1,c) - channelMin) / (channelMax - channelMin);
                if k == 1
                    imwrite(sliceNorm, strcat(plotTitle, '.tif'));    
                else
                    imwrite(sliceNorm, strcat(plotTitle, '.tif') ,'WriteMode','append');    
                end
                
                
                
            end
            
            
        end
        figure(radFig)
          
        subplot(numConditions, 1, cn);
        leg = legend(sprintfc('%d',1:numConditions));
        title(leg,'channel');
    end

 


end

writetable(colocStatsMean,[dataPath 'colocStats.csv']);