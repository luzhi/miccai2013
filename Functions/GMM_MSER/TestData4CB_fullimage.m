function [ varargout ] = TestData4CB_fullimage( gmm_model_clump, gmm_model_others, imSet, index_im2selected, imSize )
% [ gmm_post ] = TestData4NCB_fullimage( gmm_model_nuclei, gmm_model_others, imSet, index_im2selected )
%
%   Compute the posterior probabilities (nuclei, clump and background) of
%   each pixel over the whole image using the GMM model of the two guys we
%   got from the 'TrainGMM2NCB(...)'.
%
%   Input:
%       gmm_model_clump: GMM model for clump (mean, cov, prior, etc.)
%       gmm_model_others: GMM model for background (mean, cov, prior, etc.)
%
%   Output:
%       gmm_post: a cell type of posterior probabilites of the two guys
%                 of each image
%           gmm_post{i,1}(2,:) - posterior of clump
%           gmm_post{i,1}(3,:) - posterior of background
%
%       ================================================================
%       Targets: P(y = Nuclei|X), P(y = Clump|X) and P(y = Background|X)
%       ================================================================

    % Parameters
%     imSize = 1024; %512;
    MAPPING = getmapping(8,'riu2');
    
    imTest = [];
    for i = 1:length(index_im2selected)
        im = imSet{index_im2selected(i),1};
        imGradient = sobelgradient(im);
        J = lbp(im, 1, 8, MAPPING);
        imLBP = zeros(size(im));
        imLBP(2:(imSize - 1),2:(imSize - 1)) = J;
        
        % reshape the im, imGradient and imLBP from 1024x1024 to 1x1048576
        imGrayShape2Cols = reshape(im, 1, size(im,1)*size(im,2));
        imGradientShape2Cols = reshape(imGradient, 1, size(imGradient,1)*size(imGradient,2));
        imLBPShape2Cols = reshape(imLBP, 1, size(imLBP,1)*size(imLBP,2));
        
        imTest.data{i,1} = [ double( imGrayShape2Cols ) / 255; ...
            double( imGradientShape2Cols) / double( max(imGradientShape2Cols) ); ...
            double( imLBPShape2Cols / 255 ) ];        
%         imTrains.inds_nu{i,1} = 1:(size(im,1) * size(im,2));
    end
    
    gmm_post = cell(length(index_im2selected),1);
    for i = 1:length(index_im2selected)
        % Test with training data
        prob_XY_clump = pdfgauss( imTest.data{i,1}, gmm_model_clump );
        prob_XY_background = pdfgauss( imTest.data{i,1}, gmm_model_others );    

        prior_Y_clump = gmm_model_clump.Prior;
        prior_Y_background = gmm_model_others.Prior;
        
        post_Y_clump = (prob_XY_clump * prior_Y_clump) ./ ( prob_XY_clump * prior_Y_clump + prob_XY_background * prior_Y_background );
        post_Y_background = (prob_XY_background * prior_Y_background) ./ ( prob_XY_clump * prior_Y_clump + prob_XY_background * prior_Y_background );

        gmm_post{i,1} = zeros(2,imSize*imSize);
        gmm_post{i,1}(1, :) = post_Y_clump(1,:);
        gmm_post{i,1}(2, :) = post_Y_background(1,:);
    end
    
    varargout{1,1} = gmm_post;
end

