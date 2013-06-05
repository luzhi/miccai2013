function maskEllipse = drawEllipseOnCell4Expanding(nucleiMask, enlargefactor)
%
%
%   Draw ellipse around the nuclei mask
    
    maskEllipse = zeros(size(nucleiMask));
    [y,x] = find(nucleiMask == 1);
    
    X = [x,y];
    Mu = mean( X );
    XminusMu = bsxfun(@minus, X, Mu);
    
    [eigenvectorMatrix, eigenValues] = eig( XminusMu' * XminusMu ./ (size(X,1)-1) );     % covariance
    [eigenValues, order] = sort(diag(eigenValues), 'descend');
    eigenValues = 4 * diag(eigenValues); % cnstant that I cannot exaplain - experimental result (don't trust it too much)
    eigenvectorMatrix = eigenvectorMatrix(:, order);
    
    if eigenvectorMatrix(1,1) == 1 & eigenvectorMatrix(2,2) == 1
        enlargefactor = enlargefactor * 1.5;
    end
    %tt
%     enlargefactor
    
    theta = linspace(0,2*pi,100);
    eCircle = enlargefactor * [cos(theta) ; sin(theta)];        %# unit circle
    eigenvectorMatrixScale = eigenvectorMatrix * sqrt(eigenValues);               %# scale eigenvectors
    eCircle = bsxfun(@plus, eigenvectorMatrixScale*eCircle, Mu'); %#' project circle back to orig space
    
    maskEllipse = roipoly(maskEllipse, eCircle(1,:), eCircle(2,:));
%     maskEllipse = double((maskEllipse > 0).*(bwdist(maskEllipse < 0)- 0.5) - (maskEllipse < 0).*(bwdist(maskEllipse > 0)- 0.5));
%     figure; imshow(maskEllipse);
%     pause();
end