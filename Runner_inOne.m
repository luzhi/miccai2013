function Runner_inOne(dataset_name)

%clear;
%close all;
%clc;
% warning off;
%pause(0.5);

%==================================
%   Add folders/subfolders under
%   current folder
%==================================
currentfolder = pwd;
addpath(genpath(currentfolder))

%%
%==========================
% Parameters Configuration
%========================== 

%==========================
%      GMM Refinement
%==========================
wTrust = 0;
iter_in = 6;
iter_out = 1;
alfa = -2.5; %-2.5
lambda = 5;

%==========================
%    Nuclei Refinement
%==========================
wTrust_rawNuclei = 0;
iter_in_rawNuclei = 10;
iter_out_rawNuclei = 2;
alfa_rawNuclei = 5;
lambda_rawNuclei = 4;
 
%=========================
%  Nuclei Selection Rule
%=========================
maxNucleiEccentricity = 0.9;
minNucleiArea = 100;
maxNucleiClumpAreaRatio = 0.1;

%========================
%  Refine Initial Guess
%========================
iter_in_ellipse = 10;
iter_out_ellipse = 40;
alfa_ellipseSet = [-5]';
lambda_ellipseSet = [4]';
gamma_ellipseSet = [0.000001]';

Hmin =  5000; % minimal area
Hmax = 26000; % maximal area

%==========================
%        Joint LSF
%==========================
iter_in_extent = 4;
iter_out_extent = 4;
alfa_extentSet = [-5]';
lambda_extentSet = [4]'; 
gamma_extentSet = [0.15]'; 
zita_extentSet = [3]';
omega_extentSet = [4.1]';

Hmin_extent =  5000; % minimal area
Hmax_extent = 26000; % maximal area

%==========================
%    DATA FOLDERS INFO.
%==========================
%=================================
% Parameters for Training dataset
%=================================
if strcmp(dataset_name, 'Train')
		imSize = 512;
		imNum = 3;
		imFilePath = 'ims/TrainingSyntheticImages/';
    storageCommonPath = 'Train/Train_Common/';
    storageInitial = 'Train/Train_Initial/';
    storageExtent = 'Train/Train_Extent/';
end

%=================================
%   Parameters for Test dataset
%=================================
if strcmp(dataset_name, 'Test')
		imSize = 512;
		imNum = 15;
		imFilePath = 'ims/TestSyntheticImages/';
    storageCommonPath = 'Test/Test_Common/';
    storageInitial = 'Test/Test_Initial/';
    storageExtent = 'Test/Test_Extent/';
end

%=================================
%   Parameters for Real dataset
%=================================
if strcmp(dataset_name, 'EDF')
		imSize = 1024;
		imNum = 4;
		imFilePath = 'ims/EDF/';
    storageCommonPath = 'EDF/EDF_Common/';
    storageInitial = 'EDF/EDF_Initial/';
    storageExtent = 'EDF/EDF_Extent/';
end

%=================================
% Parameters for Rebuttal dataset
%=================================
if strcmp(dataset_name, 'Rebuttal')
		imSize = 512;
		imNum = 18;
		imFilePath = 'ims/Images_rebuttal/';
    storageCommonPath = 'Rebuttal/Rebuttal_Common/';
    storageInitial = 'Rebuttal/Rebuttal_Initial/';
    storageExtent = 'Rebuttal/Rebuttal_Extent/';
end

%%
%=========================
%    LOAD IMAGES DATA
%=========================
imSet = cell(imNum,1);
% imGTSet = cell(imNum,1);

%=================================
%  Read images of synthetic data
%================================= 
if strcmp(dataset_name, 'Train') || strcmp(dataset_name, 'Test') || strcmp(dataset_name, 'Rebuttal')
	 for i = 1:imNum
		   imSet{i,1} = imread(strcat(imFilePath, 'Cell', num2str(i),'.png'));
% 		   imGTSet{i,1} = imread(strcat(imFilePath, 'Cell', num2str(i),'_GTMask.png'));
	 end
end

%================================
%    Read images of real data
%================================ 
if strcmp(dataset_name, 'EDF')
	for i = 1:imNum
		  imSet{i,1} = imread(strcat(imFilePath, 'EDF00', num2str(i - 1),'.png'));
% 		  imGTSet{i,1} = imread(strcat(imFilePath, 'EDF00', num2str(i - 1),'_GTMask.png'));
	end
end

%% 
%=============================================================
% Raw Clumps boundaries by Convex Hull & Level Set refinement
%=============================================================
try
    load(strcat('Variables/', storageCommonPath, 'RawClump.mat'), 'RawClumpMaskSet');
	fprintf('Step - 1 Raw Clumps boundaries by Convex Hull & Level Set refinement...done!\n');
catch
	fprintf('Step - 1 Raw Clumps boundaries by Convex Hull & Level Set refinement...\n');
    
    RawClumpMaskSet = cell(imNum, 1);
    for i = 1:length(imSet)
				tic;
        im = imSet{i,1};
        [ LSMask ] = Preprocessing_multiLS( im );
 
        RawClumpMaskSet{i,1} = ~logical(LSMask);
        t_Convex(i) = toc;
        
        fprintf('Image %d\n', i);
    end
 		
	fprintf('done!\n');
    
    %===========================
    %       SAVE MAT FILE
    %===========================
    save(strcat('Variables/', storageCommonPath, 'RawClump.mat'), ...
        'RawClumpMaskSet', 't_Convex');
end
 
%% 
%=================================================================
% Adaptive GMM-based Training/Testing for accurate clump boundary
%=================================================================
try
    load(strcat('Variables/', storageCommonPath, 'AccurateClumpLevelSet.mat'),  'imCBMaskSet', 'imCBMaskSetRefined', 'index_im2train', 'index_im2test', 'gmm_model_clump', 'gmm_model_background', 'gmm_post_training');
    fprintf('Step - 2 GMM-based Training/Testing for accurate clump boundary...done!\n');
catch
    fprintf('Step - 2 GMM-based Training/Testing for accurate clump boundary...\n');
    
	tic;
    
    imCBMaskSet = cell(imNum,1);
    % Train the adaptive GMM model by all the images
    im_index2train = 1:imNum;
    imCBMaskSet = GenerateClumpBackgroundMask(imSet, RawClumpMaskSet, imNum);
    
    %=================================
    % Semi-supervised Learning by GMM
    %=================================
    for iter_GMM = 1:10
        if iter_GMM > 1
            imCBMaskSet = GenerateClumpBackgroundMask(imSet, imCBMaskSet, imNum);
        end
        [ gmm_model_clump, gmm_model_background ]  = TrainGMM2CB( imSet, imCBMaskSet, im_index2train, imSize );
 
        for i = 1:imNum
            im = imSet{1,1};         
            [ gmm_post_testing ] = TestData4CB_fullimage( gmm_model_clump, gmm_model_background, imSet, i, imSize );
            [ imCBMask ] = ComputeConfidenceCB(gmm_post_testing{1,1}, imSize);
            imCBMaskSet{i,1} = imCBMask;
 
        end
        
        disp(strcat('Trained ', num2str(iter_GMM)));
    end

    %=================================
    %       Refine by Level Set
    %=================================
    imCBMaskSetRefined = cell(imNum,1);
    for i = 1:imNum   
        imCBMask = imCBMaskSet{i,1};
        %+------------------------------------+
        %| Empirically Noise Removing Method. |
        %|    The usage of Level Set here     |
        %|        is different from           |
        %|       that in segmentation.        |
        %+------------------------------------+
        imCBMask = cleanNoiseRegionsByLevelSet(im, imCBMask, iter_in, iter_out, alfa, lambda, wTrust );
        imCBMask = im2bw(imCBMask);
        
        imCBMaskSetRefined{i,1} = imCBMask;
    end

    t_GMM = toc;
    fprintf('done!\n');

    %===========================
    %       SAVE MAT FILE
    %===========================
    save(strcat('Variables/', storageCommonPath, 'AccurateClumpLevelSet.mat'), ...
        'imCBMaskSet', 'imCBMaskSetRefined', 'gmm_model_clump', 'gmm_model_background', 't_GMM');
end
 
%%
%=================================
%    Using MSER to find nuclei
%================================= 
try
    load(strcat('Variables/', storageCommonPath, 'RawNucleiMask.mat'), 'imMaskSet4RawNucleiCandidatesAfterLevelSet');
	fprintf('To find nuclei...done!\n');
catch
	fprintf('To find nuclei...\n');
    
    imMaskSet4RawNucleiCandidatesAfterLevelSet = cell(length(imSet),1);
 
    %==================================
    % Using MSER to find raw candidate
    % regions for nuclei
    %==================================
    for k = 1:length(imSet)
        tic;
        
        I = imSet{k,1}; 
        
        regions = detectMSERFeatures(I, 'RegionAreaRange', [100,600]);
        rawNucleiCandidate = zeros(size(I));
        for i = 1:length(regions)    
            pixelsInRegion = regions(i,1).PixelList;
            for j = 1:size(regions(i,1).PixelList,1)
                IIx = pixelsInRegion(j,1);
                IIy = pixelsInRegion(j,2);
                rawNucleiCandidate(IIy,IIx) = 1;
            end
        end
 
        rawNucleiCandidateLogical = logical(rawNucleiCandidate);
 
        %=======================================
        %        Refinement by Level Set
        % (remove noises around nuclei regions)
        %=======================================
        imMask4RawNucleiCandidatesAfterLevelSet = cleanNoiseRegionsByLevelSet(I, rawNucleiCandidateLogical, ... 
                                                                                 iter_in_rawNuclei,...
                                                                                 iter_out_rawNuclei, alfa_rawNuclei,...
                                                                                 lambda_rawNuclei, wTrust_rawNuclei );
 
        imMask4RawNucleiCandidatesAfterLevelSet = im2bw(imMask4RawNucleiCandidatesAfterLevelSet); 
 
        imMaskSet4RawNucleiCandidatesAfterLevelSet{k,1} = imMask4RawNucleiCandidatesAfterLevelSet;
        
        t_raw_nuclei_mask(k) = toc;
    end
 
    fprintf('done!\n');
    
    %===========================
    %       SAVE MAT FILE
    %===========================
    save(strcat('Variables/', storageCommonPath, 'RawNucleiMask.mat'), ... 
        'imMaskSet4RawNucleiCandidatesAfterLevelSet', 't_raw_nuclei_mask');
end
 
%%
%=============================
% Using Rules to find nucleus
% (from RawNucleiCandidates)
%============================= 
try
    load(strcat('Variables/', storageCommonPath, 'NucleiMask.mat'), 'nucleiMaskSet');
    fprintf('Rule-based Nuclei Selection...done!\n');
catch
    fprintf('Rule-based Nuclei Selection...\n');
    
    %======================================
    % Use rules to remove false nucleus...
    %======================================
    nucleiMaskSet = cell(imNum,1);
    for i = 1:imNum
        tic;
        
        clumpBackgroundMask = imCBMaskSetRefined{i,1};
        nucleicandidatesMask = ~imMaskSet4RawNucleiCandidatesAfterLevelSet{i,1};
        nucleiMask = zeros(size(nucleicandidatesMask));
        
        regionNuleiStats = regionprops(nucleicandidatesMask, 'Eccentricity', 'PixelIdxList', 'Area');
        
        for k = 1:length(regionNuleiStats)
            
            areaRatio = computeAreaRatio4NucleiClump(regionNuleiStats(k,1).Area, ...
                        regionNuleiStats(k,1).PixelIdxList, clumpBackgroundMask);
            if regionNuleiStats(k,1).Eccentricity < maxNucleiEccentricity && ...
                    regionNuleiStats(k,1).Area > minNucleiArea &&...
                    areaRatio < maxNucleiClumpAreaRatio
                nucleiMask( regionNuleiStats(k,1).PixelIdxList(:) ) = 1;
            end
        end
        
        nucleiMaskSet{i,1} = nucleiMask;
        
        t_refined_nuclei_mask(i) = toc;
    end

    fprintf('done!\n');
    
    %===========================
    %       SAVE MAT FILE
    %===========================
    save(strcat('Variables/', storageCommonPath, 'NucleiMask.mat'), ... 
        'nucleiMaskSet', 't_refined_nuclei_mask');
end
 
%%
%==============================================
%    Use Level Set-based method to find out 
%   the individual cell boundary from a clump
%==============================================
try
    load('xx');
catch
    fprintf('Initial LSF evolution...\n');
    
    for alphaID = 1:length(alfa_ellipseSet)				
        for lambdaID = 1:length(lambda_ellipseSet)						
            for gammaID = 1:length(gamma_ellipseSet)
                
                alfa_ellipse = alfa_ellipseSet(alphaID,1);
                lambda_ellipse = lambda_ellipseSet(lambdaID,1);
                gamma_ellipse = gamma_ellipseSet(gammaID,1);

                phi_masks_set = cell(imNum,1);

                for i = 1:imNum
                    tic;                    
                    fprintf('Image %d.png\n', i);
                    
                    im = imSet{i,1};
                    clumpRefinedMask = imCBMaskSetRefined{i,1};
                    nucleiMask = logical(nucleiMaskSet{i,1});

                    % Get clump & nuclei regions pixels' index
                    clumpStats = regionprops(clumpRefinedMask, 'PixelIdxList');
                    nucleiStats = regionprops(nucleiMask, 'PixelIdxList');


                    cellMasks = cell(1,1);
                    nucleiCounter = 1;

                    for k = 1:length(clumpStats)
                        phi_nucleus_inside_clump = zeros(size(im));
                        clumpMask = zeros(size(im));
                        clumpMask( clumpStats(k,1).PixelIdxList ) = 1;
                        clumpMask = logical(clumpMask);
                        I = im;
                        I( clumpMask ~= 1 ) = 0;


                        % find the nucleus inside this clump
                        insideNucleiNum = 0;
                        for n = 1:length(nucleiStats)                
                            isNucleiInsideClump = isequal(intersect(nucleiStats(n,1).PixelIdxList(:),...
                                                  clumpStats(k,1).PixelIdxList(:)),...
                                                  nucleiStats(n,1).PixelIdxList(:));
                            if isNucleiInsideClump == 1
                                phi_nucleus_inside_clump( nucleiStats(n,1).PixelIdxList ) = 1;
                                insideNucleiNum = insideNucleiNum + 1;
                            end
                        end

                        % if the number of nuclei in this clump is "insideNucleiNum == 1"
                        % store its clump curve as the cytoplasm curve
                        if insideNucleiNum == 1
                            cellMasks{nucleiCounter,1} = logical(clumpMask);
                            nucleiCounter = nucleiCounter + 1;
                            continue;
                        end
                        
                        if insideNucleiNum == 0
                            fprintf('\tFound 0 nuclei in clump %d of Image %d\n', k, i);
                            continue;
                        end

                        % otherwise, evolve the curve by the nucleus inside the clump
                        phi_nucleus_inside_clump = logical(phi_nucleus_inside_clump);        

                        % use ellipse term & Level Set to find out the boundary of
                        % individual cell inside this clump one by one
                        phiStats = regionprops(phi_nucleus_inside_clump, 'PixelIdxList');

                        % tt
                        fprintf('\tFound %d nucleus in clump %d of Image %d\n', length(phiStats), k, i);

                        % detect individual cell boundary one nuclei by another
                        % (NOT simultaneously)            
                        for m = 1:length(phiStats)
                            phi = zeros(size(im));
                            phi( phiStats(m,1).PixelIdxList ) = 1;
                            phi = logical(phi);

                            phi = drawEllipseOnCell(phi, 1);                         
                            
                            
                            %=======================
                            %     LSF evolution
                            %=======================														
                            phi = ellipseLevelSet(I, clumpMask, phi, iter_in_ellipse,...
                                                  iter_out_ellipse, alfa_ellipse,...
                                                  lambda_ellipse, gamma_ellipse, Hmin, Hmax);

                            cellMasks{nucleiCounter,1} = phi;
                            nucleiCounter = nucleiCounter + 1;
                        end
                    end
                    phi_masks_set{i,1} = cellMasks;
                    t_Initial_LSF(i) = toc;
                end
                
                %===========================
                %       SAVE MAT FILE
                %===========================       
                save(strcat('Variables/',storageInitial, '/Initial_LSF_alfa', num2str(alfa_ellipse),'_lambda', num2str(lambda_ellipse),...
                            '_gamma', num2str(gamma_ellipse), '_Hmin', num2str(Hmin),...
                            '_Hmax', num2str(Hmax), '_iterIn', num2str(iter_in_ellipse), '_iterOut', num2str(iter_out_ellipse), '.mat'),...
                            'phi_masks_set', 't_Initial_LSF');
            end
        end
    end
    fprintf('done!\n');
end


%%
%==============================
%           Joint LSF
%==============================
try
    load('xxxx');
catch
    for alfa_extentID = 1:length(alfa_extentSet)             
        for lambda_extentID = 1:length(lambda_extentSet)           
            for gamma_extentID = 1:length(gamma_extentSet)                  
                for zita_extentID = 1:length(zita_extentSet)                
                    for omega_extentID = 1:length(omega_extentSet)
                        
                        alfa_extent = alfa_extentSet(alfa_extentID);
                        lambda_extent = lambda_extentSet(lambda_extentID);
                        gamma_extent = gamma_extentSet(gamma_extentID);
                        zita_extent = zita_extentSet(zita_extentID);
                        omega_extent = omega_extentSet(omega_extentID);
                        
                        fprintf('Alpha = %f, Lambda = %f, Gamma = %f, Zita = %f, Omega = %f\n', ...
								alfa_extent, lambda_extent, gamma_extent, zita_extent, omega_extent);

                        phi_refinedCytoplasms_masks_set = cell(imNum,1);   % store the refined individual cells

                        for i = 1:imNum
                            tic;

                            im = imSet{i,1};
                            rawCytoplasmsPhiSet = phi_masks_set{i,1};  % clumps of individual cells in previous step
                            % covert double phi to logical phi
                            for k = 1:length(rawCytoplasmsPhiSet)
                                if islogical(rawCytoplasmsPhiSet{k,1}) == 0
                                    rawCytoplasmsPhiSet{k,1} = ~im2bw(rawCytoplasmsPhiSet{k,1});
                                end
                            end

                            clumpsMask = imCBMaskSetRefined{i,1};   % clumps background mask
                            clumpsStats = regionprops(clumpsMask, 'PixelIdxList');  % pixels lists of each clump
                            phiCytoplasmsMasksInOneClump = cell(length(rawCytoplasmsPhiSet),1);  %s store the refined cells' phis of a clump

                            for k = 1:length(rawCytoplasmsPhiSet)
                                phi_1 = rawCytoplasmsPhiSet{k,1};
                                phi_1_Idx = find(phi_1 == 1);
                                clumpContainPhi1 = zeros(size(im));
                                neighborsOfPhi_1 = cell(1,1);
                                neighborNum = 1;
                                I = im;

                                % find the clump contains this cytoplasm
                                for j = 1:length(clumpsStats)
                                    clumpPixelsIdx = clumpsStats(j,1).PixelIdxList;
                                    isNotInteract = isempty(intersect(phi_1_Idx, clumpPixelsIdx));

                                    if isNotInteract == 0
                                        clumpContainPhi1(clumpPixelsIdx) = 1;
                                        clumpContainPhi1 = logical(clumpContainPhi1);
                                        I( clumpContainPhi1 ~= 1 ) = 0;
                                    end
                                end

                                % find the neighbour cytoplasms of this cytoplasm
                                for j = 1:length(rawCytoplasmsPhiSet)
                                    if j == k
                                        continue;
                                    end
                                    phi2TempImg = rawCytoplasmsPhiSet{j,1};
                                    phi2PixelsIdx = find(phi2TempImg == 1);

                                    isNotInteract = isempty(intersect(phi_1_Idx, phi2PixelsIdx));

                                    if isNotInteract == 0
                                        neighborsOfPhi_1{neighborNum,1} = phi2TempImg;
                                        neighborNum = neighborNum + 1;
                                    end
                                end

                                % for single-cell, push it into the clumps stack
                                if length(neighborsOfPhi_1) == 1 & isempty(neighborsOfPhi_1{1,1}) == 1
                                    phiCytoplasmsMasksInOneClump{k,1} = phi_1;
                                    continue;
                                end

                                % run level set for each neighbour cytoplasm of this cytoplasm            
                                phi_extent = phi_1;
                                for j = 1:length(neighborsOfPhi_1)
                                    phi_2 = neighborsOfPhi_1{j,1};
                                    if islogical(phi_extent) == 0
                                        phi_extent = ~im2bw(phi_extent);
                                    end

                                    %======================
                                    %    LSF evolution
                                    %======================
                                    phi_extent = overlapExtentLevelSet(im, clumpContainPhi1, phi_extent, phi_2, iter_in_extent,...
                                                                       iter_out_extent, alfa_extent,...
                                                                       lambda_extent, gamma_extent, zita_extent, omega_extent, Hmin_extent, Hmax_extent);
                                end

                                phiCytoplasmsMasksInOneClump{k,1} = phi_extent;
                            end
                            phi_refinedCytoplasms_masks_set{i,1} = phiCytoplasmsMasksInOneClump;
                                    t_final_LSF(i) = toc;
                        end

                        fprintf('done!\n');
                        
                        %===========================
                        %       SAVE MAT FILE
                        %===========================
                        save(strcat('Variables/',storageExtent,'/Final_LSF_alfa', num2str(alfa_extent),'_lambda', num2str(lambda_extent),...
                                        '_gamma', num2str(gamma_extent), '_zita', num2str(zita_extent), '_omega', num2str(omega_extent), '_Hmin', num2str(Hmin_extent),...
                                        '_Hmax', num2str(Hmax_extent), '_iterIn', num2str(iter_in_extent), '_iterOut', num2str(iter_out_extent), '.mat'),...
                                        'phi_refinedCytoplasms_masks_set', 't_final_LSF');
                    end
                end
            end
        end
    end
end
