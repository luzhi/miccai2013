function [ LSMask ] = Preprocessing_multiLS( im)
% Preprocessing - ...
% X,Y,K and V for segmentation over each convex hull
% K - indicator of the hulls
% V - the volume (number of pixels) in each hull
% tic;
    %% parameters
    handles.ad_k = 0.05;
    handles.ad_iter = 5;
    handles.conhull_thr = 400;
    handles.thr_tinyfrags = 400;
    
    %% common steps
    % denoise
    img_denoised = anisodiff(im, handles.ad_iter, ...
                                    handles.ad_k, 0.25, 1);
    img_denoised = uint8(img_denoised);
    % quick shift - superpixels
%     disp('Meanshift...');
%     disp(' ');
    [imSegbySP,SPLabels] = vl_quickseg(img_denoised,0.5,2,10);
    
    %
    imGrdbyMorph = morphGrad(imSegbySP);
    bwImGrdbyMorph = bwGrdByThr(imGrdbyMorph,0.03);

    % skeleton the lines
    bwImGrdbyMorph = bwmorph(bwImGrdbyMorph,'thin');  
        
    [meanShiftLabels, num] = bwlabel(bwImGrdbyMorph, 8);
    % X,Y,K and V for segmentation over each convex hull
    % K - indicator of the hulls
    % V - the volume (number of pixels) in each hull

    phi_convexHull = zeros(size(im));
    for i=1:max(max(meanShiftLabels))
        [Y,X] = find(meanShiftLabels==i);
        if (length(Y) > handles.conhull_thr)            
            [K, V]=convhull(X,Y);                        
            subImageMask = roipoly(im,X(K),Y(K));
            phi_convexHull( subImageMask == 1 ) = 1;
        end;
    end;
    
%     %tt
%     figure; imshow(im); hold on; contour(phi_convexHull, [0,0], 'r');
    
    LSMask = im2bw(phi_convexHull);
%     phi = cleanNoiseRegionsByLevelSet(im, phi_convexHull, 5, 40, 1.5, 5, 0);
%   
%     phiBWLevel = graythresh(phi);
%     LSMask = im2bw(phi,phiBWLevel);
    
%     % tt
%     figure; imshow(im); hold on; contour( LSMask, [0,0], 'b');                     
% toc;
end

