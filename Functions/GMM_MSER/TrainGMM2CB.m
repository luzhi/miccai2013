function [ varargout ] = TrainGMM2CB( imSet, imNCBMasks, index_im2train, imSize )
% [ gmm_model_nuclei, gmm_model_clump, gmm_model_background ] = TrainGMMForeBackground_noNuclei4Training( imSet, imGTSet, imClumpEdgesSet, index_im2train )
%   Train the parameters (MuK, SigmaK, PiK) of a GMM for the clump
%   (including nuclei) and the background (named as 'others')
%   
%   Features: gray + gradient + lbp
%
%   Input:
%       imSet: a set of images for training
%       imNCBMasks: a cell type of the masks of nuclei, clump and
%                   background
%                   nuclei: 0
%                   clump:  100
%                   background: 255
%       index_im2train: radomly selected images for training
%
%   Output:
%       gmm_model_clump: GMM model for clump (mean, cov, prior, etc.)
%       gmm_model_others: GMM model for background (mean, cov, prior, etc.)
%
%       ===============================================
%       Training targets: P(X|theta_C) and P(X|theta_B)
%       ===============================================
%
%   Example:
%       imSet = cell(4,1);
%       imSet{1,1} = imread('ims\EDF000.png');
%       imSet{2,1} = imread('ims\EDF001.png');
%       imSet{3,1} = imread('ims\EDF002.png');
%       imSet{4,1} = imread('ims\EDF003.png');
%       imGTSet = cell(4,1);
%       imGTSet{1,1} = imread('ims\EDF000_GTMask.png');
%       imGTSet{2,1} = imread('ims\EDF001_GTMask.png');
%       imGTSet{3,1} = imread('ims\EDF002_GTMask.png');
%       imGTSet{4,1} = imread('ims\EDF003_GTMask.png');
%
%       [ gmm_model_clump, gmm_model_others ] = ...
%       TrainGMMForeBackground_noNuclei4Training( imSet, imGTSet, imClumpEdgesSet, index_im2train );

    % Parameters
%     imSize = 1024; %512;
    MAPPING = getmapping(8,'riu2');

    tmpData4train = [];
    
    tmpData4train.clump_gray = -1;
    tmpData4train.clump_gradient = -1;
    tmpData4train.clump_lbp = -1;
    
    tmpData4train.others_gray = -1;
    tmpData4train.others_gradient = -1;
    tmpData4train.others_lbp = -1;
    
    data4train = [];    
    
    for i = 1:length(index_im2train)
        % features extraction
        im = imSet{index_im2train(i),1};
        imNCB = imNCBMasks{index_im2train(i),1};
        imGradient = sobelgradient(im);        
        J = lbp(im, 1, 8, MAPPING);
        imLBP = zeros(size(im));
        imLBP(2:(imSize - 1),2:(imSize - 1)) = J;
        
        % ROI of the image    
        
        inds_clump_train = imNCB == 100;
        
        inds_background_train = imNCB == 255;
        
        tmpData4train.clump_gray =[tmpData4train.clump_gray; im( inds_clump_train )];
        tmpData4train.clump_gradient = [tmpData4train.clump_gradient; imGradient( inds_clump_train )];
        tmpData4train.clump_lbp = [tmpData4train.clump_lbp; imLBP( inds_clump_train )];
        
        tmpData4train.others_gray = [tmpData4train.others_gray; im( inds_background_train )];
        tmpData4train.others_gradient = [tmpData4train.others_gradient; imGradient( inds_background_train )];
        tmpData4train.others_lbp = [tmpData4train.others_lbp; imLBP( inds_background_train )];
    end
    
    % the Training data set! + normalization
    
    data4train.clump = [tmpData4train.clump_gray(2:length(tmpData4train.clump_gray)), tmpData4train.clump_gradient(2:length(tmpData4train.clump_gradient)), tmpData4train.clump_lbp(2:length(tmpData4train.clump_lbp))];
    data4train.clump = double(data4train.clump);
    data4train.clump(:,1) = double( data4train.clump(:,1) / 255 );
    data4train.clump(:,2) = double( data4train.clump(:,2) / max(data4train.clump(:,2)) );
    data4train.clump(:,3) = double( data4train.clump(:,3) / 255 );
    
    data4train.clump_labels = ones(length(tmpData4train.clump_gray) - 1,1);
    
    data4train.others = [tmpData4train.others_gray(2:length(tmpData4train.others_gray)), tmpData4train.others_gradient(2:length(tmpData4train.others_gradient)), tmpData4train.others_lbp(2:length(tmpData4train.others_lbp))];
    data4train.others = double(data4train.others);
    data4train.others(:,1) = double( data4train.others(:,1) / 255 );
    data4train.others(:,2) = double( data4train.others(:,2) / max(data4train.others(:,2)) );
    data4train.others(:,3) = double( data4train.others(:,3) / 255 );
    
    data4train.others_labels = 2 * ones(length(tmpData4train.others_gray) - 1,1);

    DSTrain_clump = [];
    DSTrain_clump.X = data4train.clump(:,:)';    
    DSTrain_clump.y = data4train.clump_labels(:,:)'; 
    DSTrain_clump.name = 'Training set - Clump';
    DSTrain_clump.dim = size(data4train.clump,2);
    DSTrain_clump.num_data = size(data4train.clump,1);
    
    DSTrain_others = [];
    DSTrain_others.X = data4train.others(:,:)';
    DSTrain_others.y = data4train.others_labels(:,:)';
    DSTrain_others.name = 'Training set - Background';
    DSTrain_others.dim = size(data4train.others,2);
    DSTrain_others.num_data = size(data4train.others,1);
    
    % Estimate the paramters of GMM
    gmm_model_clump = mlcgmm(DSTrain_clump);
    gmm_model_others = mlcgmm(DSTrain_others);
    
    % Ouput estimated GMM parameters
    varargout{1,1} = gmm_model_clump;
    varargout{2,1} = gmm_model_others;
end

