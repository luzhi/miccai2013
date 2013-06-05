function [ varargout ] = ellipseLevelSet( im, clumpMask, nuclei_MASK, iter_inner, iter_outer, alfa, lambda, gamma, Hmin, Hmax )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
    
    % storage of phis for different iter_out
%     Phis = cell(1,1);
%     SaveNum = 1;
    
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
    initialLSF(nuclei_MASK == 1)=-c0;  
    phi=initialLSF;
 
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
%         gamma
        if gamma >= 0.03
            gamma = 0.03;
        end
%         if gamma >= 0.000002
%             gamma = 0.000002
%         end
        
        
        disp(strcat('             Ellipse Level Set - iteration #: ', num2str(n)));
        phi = ellipse_drlse_edge(im, phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, gamma, potentialFunction, Hmin, Hmax);
        phi_tmp = ~im2bw(phi);
        phi_1 = phi_tmp & clumpMask;
 
        phi_1 = double(phi_1);
        phi_1( phi_1 == 1 ) = -c0;
        phi_1( phi_1 == 0 ) = c0;
        phi = double((phi_1 > 0).*(bwdist(phi_1 < 0)- 0.5) - (phi_1 < 0).*(bwdist(phi_1 > 0)- 0.5));
        
        % update some parameters...        
%         gamma = gamma * 1.05;
        gamma = gamma * 1.4;
%         zita = zita * 1.1;
%  
%         if mod(n,1)==0
%             figure(99);clf;
%             imshow(im); hold on;  contour(im2bw(phi), [0,0], 'r'); hold on;
%             pause(0.5);
%         end
        
%         if mod(n,10) == 0
%             Phis{SaveNum, 1} = phi;
%             SaveNum = SaveNum + 1;
%         end
    end
    
    disp(' ');
 
    varargout{1,1} = phi;

end
