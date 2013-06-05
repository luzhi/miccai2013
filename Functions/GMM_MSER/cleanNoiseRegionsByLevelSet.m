function [ varargout ] = cleanNoiseRegionsByLevelSet( im, im_conf_nuclei_MASK, iter_inner, iter_outer, alfa, lambda, wTrust )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
    
    Img=double(im(:,:,1));
    % parameter setting
    timestep=5;  % time step
    mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
%     iter_inner=5;
%     iter_outer=40;
%     lambda=5; % coefficient of the weighted length term L(phi)
%     alfa=-1.5;  % coefficient of the weighted area term A(phi)
    epsilon=1.5; % papramater that specifies the width of the DiracDelta function

    sigma=1.5;     % scale parameter in Gaussian kernel
    G=fspecial('gaussian',15,sigma);
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.

    % initialize LSF as binary step function
    c0=2;
    initialLSF=c0*ones(size(Img));
    % generate the initial region R0 as a rectangle
    initialLSF(im_conf_nuclei_MASK == 1)=-c0;  
    phi=initialLSF;

%     figure(1);
%     mesh(-phi);   % for a better view, the LSF is displayed upside down
%     hold on;  contour(phi, [0,0], 'r','LineWidth',2);
%     title('Initial level set function');
%     % VIEW([-80 35]);
% 
%     figure(2);
%     imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%     title('Initial zero level contour');
%     pause(0.5);

    potential=2;  
    if potential ==1
        potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
    elseif potential == 2
        potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
    else
        potentialFunction = 'double-well';  % default choice of potential function
    end


    % start level set evolution
    for n=1:iter_outer
        disp(strcat('    Level Set - iteration #: ', num2str(n)));
        phi = drlse_edge_for_clump_background(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction, wTrust);
%         if mod(n,2)==0
%             figure(2);
%             imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%         end
    end

    % refine the zero level contour by further level set evolution with alfa=0
    disp(strcat('    Level Set - iteration # (refinement loop): ', num2str(iter_outer + 1)));
    disp(' ');
    alfa=0;
    iter_refine = 10;
    phi = drlse_edge_for_clump_background(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction, wTrust);

%     finalLSF=phi;
%     figure(2);
%     imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%     hold on;  contour(phi, [0,0], 'r');
%     str=['Final zero level contour, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
%     title(str);
    varargout{1,1} = phi;
end

