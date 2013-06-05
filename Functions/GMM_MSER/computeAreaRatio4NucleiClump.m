function [ NucleiClumpAreaRatio ] = computeAreaRatio4NucleiClump(nucleiRegionArea, nucleiRegionPixelIdxList, imClumpMask)
%
%
%   Compute a set of area ratios between nuclei regions and its clump
%   regions
    
    NucleiClumpAreaRatio = -1;
    clumpStats = regionprops(imClumpMask, 'Area', 'PixelIdxList');


    for k = 1:length(clumpStats)
        isNucleiBelongToClump = isequal(intersect(nucleiRegionPixelIdxList(:), clumpStats(k,1).PixelIdxList(:)), nucleiRegionPixelIdxList(:));
        if isNucleiBelongToClump == 1
            NucleiClumpAreaRatio = double( nucleiRegionArea / clumpStats(k,1).Area );
            break;
        end
    end

    if NucleiClumpAreaRatio == -1;
        NucleiClumpAreaRatio = 0.000001;
    end
end