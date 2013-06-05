function [ varargout ] = ComputeConfidenceCB( gmm_post, imSize )
% [ varargout ] = ComputeConfidenceNCB( gmm_post )
%   Depending on the posterior obtained by GMM, calculate the labels for
%   each pixels by confidence
%
%   Method: compute ratio between any pair of NCB probabilites
%
%   Input:
%       gmm_post: elements of a cell of posteriors of training images (NCB)
%
%   Output:
%       varargout{2,1} - confidence between clump and background

%     imSize = 1024; %512;
    im_clump = reshape(gmm_post(1,:), imSize,imSize);
    im_background = reshape(gmm_post(2,:), imSize, imSize);
    
    % confidence between clump and background
    im_conf_clump_background = exp( im_clump - im_background );
    im_conf_clump_background = im_conf_clump_background / max(max( im_conf_clump_background ));
    im_conf_clump_background = uint8( im_conf_clump_background * 255 );
    im_conf_clump_background = im2bw( im_conf_clump_background );
    
    % confidence between nuclei and others (clump + background)
    % [ in progress... ]
    
   varargout{1,1} = im_conf_clump_background;
end

