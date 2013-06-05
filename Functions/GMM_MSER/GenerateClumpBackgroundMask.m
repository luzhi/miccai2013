function [ imClumpBackgroundMaskSet ] = GenerateClumpBackgroundMask( imSet, imGTSet, RawClumpMaskSet, imNum )
% Generate the mask containing nuclei, clump and background using Andrew's
% annotated nuclei mask and our previous generated clump
%   Detailed explanation goes here
    
    imClumpBackgroundMaskSet = cell(imNum,1);
    for i = 1:imNum
        im = imSet{i,1};
        % draw background
        imClumpBackground = 255 * ones(size(im));

        % draw clump
        ClumpMask4SingleImage = RawClumpMaskSet{i,1};
        imClumpBackground(sub2ind(size(imClumpBackground), find(ClumpMask4SingleImage == 1))) = 100;

        imClumpBackgroundMaskSet{i,1} = uint8(imClumpBackground);
    end
end

