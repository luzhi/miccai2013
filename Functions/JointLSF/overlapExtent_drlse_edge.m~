function phi = overlapExtent_drlse_edge(Img, phi_0, phi_mask_another, g, lambda, mu, alfa, epsilon, timestep, iter, gamma, zita, omega, potentialFunction, Hmin, Hmax)
%  This Matlab code implements an edge-based active contour model as an
%  application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper:
%
%      C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation", 
%        IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
%
%  Input:
%      phi_0: level set function to be updated by level set evolution
%      g: edge indicator function
%      mu: weight of distance regularization term
%      timestep: time step
%      lambda: weight of the weighted length term
%      alfa:   weight of the weighted area term
%      epsilon: width of Dirac Delta function
%      iter: number of iterations
%      potentialFunction: choice of potential function in distance regularization term. 
%              As mentioned in the above paper, two choices are provided: potentialFunction='single-well' or
%              potentialFunction='double-well', which correspond to the potential functions p1 (single-well) 
%              and p2 (double-well), respectively.%
%  Output:
%      phi: updated level set function after level set evolution
%
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.engr.uconn.edu/~cmli/

phi=phi_0;
phi_mask_another = ~phi_mask_another;


% Hmin = 9000;
% Hmax = 15000;
% iRecorder = [];
        
[vx, vy]=gradient(g);
gray = Img;

for k=1:iter
    phi=NeumannBoundCond(phi);
    [phi_x,phi_y]=gradient(phi);
    s=sqrt(phi_x.^2 + phi_y.^2);
    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
    Ny=phi_y./(s+smallNumber);
    curvature=div(Nx,Ny);
    if strcmp(potentialFunction,'single-well')
        distRegTerm = 4*del2(phi)-curvature;  % compute distance regularization term in equation (13) with the single-well potential p1.
    elseif strcmp(potentialFunction,'double-well');
        distRegTerm=distReg_p2(phi);  % compute the distance regularization term in eqaution (13) with the double-well potential p2.
    else
        disp('Error: Wrong choice of potential function. Please input the string "single-well" or "double-well" in the drlse_edge function.');
    end        
    
    
    diracPhi=Dirac(phi,epsilon);
    phi1_Size = phiArea(phi);
    phi2_Size = phiArea(phi_mask_another);
   
%   Area Term
    if phi1_Size <= Hmax && phi1_Size >= Hmin
        areaTerm = 0;
    else if phi1_Size < Hmin
            areaTerm= -1 * (g * phi1_Size - Hmin) .* diracPhi.*g; % balloon/pressure force
        else if phi1_Size > Hmax
            areaTerm= 1 * (g * phi1_Size - Hmax) .* diracPhi.*g; % balloon/pressure force
            end
        end
    end
    
%   Ellipse Shape Term
    phi_ellipse = drawEllipseOnCell4Expanding(~im2bw(phi), 1);    
    phi_ellipse = double(phi_ellipse);
%     % try 1
%     if k == 1
%         phi_ellipse( phi_ellipse == 1 ) = -2;
%         phi_ellipse( phi_ellipse == 0 ) = 2;
%     else
%         phi_ellipse( phi_ellipse == 1 ) = -2;
%         phi_ellipse( phi_ellipse == 0 ) = 2;
%         phi_ellipse = double((phi_ellipse > 0).*(bwdist(phi_ellipse < 0)- 0.5) - (phi_ellipse < 0).*(bwdist(phi_ellipse > 0)- 0.5)); %??
%     end

    % try 2
    phi_ellipse( phi_ellipse == 1 ) = -2;
    phi_ellipse( phi_ellipse == 0 ) = 2;
    phi_ellipse = double((phi_ellipse > 0).*(bwdist(phi_ellipse < 0)- 0.5) - (phi_ellipse < 0).*(bwdist(phi_ellipse > 0)- 0.5));
    
    
    diracEllipsePhi = Dirac(phi_ellipse, epsilon);
    ellipseTerm = diracEllipsePhi .* g;
    
    
%   Edge Term
    edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
    
%   Overlapping Area Extent Term
%     phi_mask_another = ~phi_mask_another;
    H_Phi1 = Heaviside(phi);
    H_Phi2 = Heaviside(phi_mask_another);
    IntersectH_Phi = H_Phi1 .* H_Phi2;
    IntersectArea = sum(IntersectH_Phi(:) == 1);
    
    ratioIntersectArea = IntersectArea / phi1_Size;
    
    if ratioIntersectArea >= 0.10
        phi_overlap = -diracPhi .* g .* ( H_Phi2 * phiArea(g.*phi) - phiArea(g.*IntersectH_Phi)) ./ ( phiArea(g.*phi)^2);
%         phi_overlap = diracPhi .* g .* ( H_Phi2 * phiArea(g.*phi) - phiArea(g.*IntersectH_Phi)) ./ ( phiArea(g.*phi)^2);
    else
         phi_overlap = 0;
%        phi_overlap = diracPhi .* g .* ( H_Phi2 * phiArea(g.*phi) - phiArea(g.*IntersectH_Phi)) ./ ( phiArea(g.*phi)^2);
    end
    % tt
%     disp(strcat('phi1_size = ', num2str(phi1_Size), '    InteractArea = ', num2str(IntersectArea), '    ratio = ', num2str(ratioIntersectArea)));

    % gray term
%     y_grayTerm = gray .* diracPhi .* (H_Phi2 * g * phi2_Size - g * IntersectArea) ./ (( g .* phi1_Size) * ( g .* phi1_Size));
    
    phi1GrayH = gray .* H_Phi1;
    IntersectGrayH = gray .* IntersectH_Phi;
    
    phi1Gray = mean(phi1GrayH(phi1GrayH ~= 0));
    IntersectGray = mean(IntersectGrayH(IntersectGrayH ~= 0));
    
    if IntersectGray <= (phi1Gray + 10)
        %grayTerm = -1 * ((gray .* diracPhi .* g .* phi1_Size - gray .* phi1_Size .* g .* diracPhi) ./ ((g .* phi1_Size) .* (g .* phi1_Size)))...
        %       + ((gray .* diracPhi .* H_Phi2 .* g .* IntersectArea - gray .* IntersectArea .* g .* diracPhi .* H_Phi2) ./ ((g .* IntersectArea) .* (g .* IntersectArea)));
        grayTerm = 0;   
        disp('darker!');
    else
        grayTerm = -(...
            - ((gray .* diracPhi .* (1-H_Phi2) .* phiArea(g.*phi.*(1-H_Phi2)) - phiArea(gray.*phi.*(1-H_Phi2)) .* g .* diracPhi .* (1-H_Phi2)) ./ (phiArea(g.*phi.*(1-H_Phi2))^2))...
               + ((gray .* diracPhi .* H_Phi2 .* phiArea(g.*IntersectH_Phi) - gray .* diracPhi .* H_Phi2 * phiArea(gray.*IntersectH_Phi)) ./...
               (phiArea(g.*IntersectH_Phi)^2)));
        disp('lighter!');
    end
    disp(strcat('phi1_size = ', num2str(phi1_Size), '    InteractArea = ', num2str(IntersectArea), '    ratio = ', num2str(ratioIntersectArea),...
        '    phi1Gray = ', num2str(phi1Gray), '    IntersectGray = ', num2str(IntersectGray)));
    
    
    if isnan(grayTerm(1))
        disp('gotcha');
    end;
    
%     phi = phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm + gamma * phi_ellipse + 10 * phi_overlap);
%     phi = phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm + gamma * phi_ellipse + zita * phi_overlap + omega * grayTerm);
    phi = phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm + gamma * ellipseTerm + zita * phi_overlap + omega * grayTerm);

%     phi = phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm + gamma * phi_ellipse);

%     phi = phi + timestep*(mu*distRegTerm + lambda*edgeTerm + 0*areaTerm + 0 * phi_ellipse + zita * phi_overlap + omega * grayTerm);

%     phi = double((phi > 0).*(bwdist(phi < 0)- 0.5) - (phi < 0).*(bwdist(phi > 0)- 0.5));  
    
    %tt
    if 0
    figure(99);
    imshow(Img,[0, 255]); hold on;  contour(phi, [0,0], 'r'); hold on;
    contour(phi_ellipse, [0,0], 'g'); % hold on; contour(phi_mask_another, [0,0], 'b');
    hold on; contour(phi_mask_another, [0,0], 'b');
    pause(0.1);
    end;
    
%     if 1
%         figure(299);
%         mesh(-phi);
%         pause(0.1);
%     end
    
%     figure(11);
%     mesh(-phi);   % for a better view, the LSF is displayed upside down
%     hold on;  contour(phi, [0,0], 'r','LineWidth',2);
%     pause(0.1);
end



function f = Heaviside(phi)
% covert distance function phi to a binary Heaviside mask
f = ~im2bw(phi);

function f = phiArea(phi)
% count the number of pixels inside phi
phiBw = ~im2bw(phi);
f = length(find(phiBw == 1));
% phiStats = regionprops(phiBw, 'Area', 'PixelIdxList');
% phiLabel = bwlabel(phiBw, 8);
% phiIdxInLabel = find(phiLabel == 1);
% for i = 1:length(phiStats)
%     isPhi = isequal(phiStats(i,1).PixelIdxList(:), phiIdxInLabel);
%     if isPhi == 1
%         f = phiStats(i,1).Area;
%         break;
%     end
% end


function f = distReg_p2(phi)
% compute the distance regularization term with the double-well potential p2 in eqaution (16)
[phi_x,phi_y]=gradient(phi);
s=sqrt(phi_x.^2 + phi_y.^2);
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
f = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi);  

function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  
