function phi = ellipse_drlse_edge(Img, phi_0, g, lambda, mu, alfa, epsilon, timestep, iter, gamma, potentialFunction, Hmin, Hmax)
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
        
[vx, vy]=gradient(g);

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
    phiSize = phiArea(phi);

    %=======================
    %       Area term
    %=======================
    if phiSize <= Hmax & phiSize >= Hmin
        areaTerm = 0;
    else if phiSize < Hmin
            areaTerm= -1 * (g * phiSize - Hmin) .* diracPhi.*g; % balloon/pressure force
        else if phiSize > Hmax
            areaTerm= 1 * (g * phiSize - Hmax) .* diracPhi.*g; % balloon/pressure force
            end
        end
    end

    %=======================
    %     Ellipse term
    %=======================
    phi_ellipse = drawEllipseOnCell4Expanding(~im2bw(phi), 1);    
    phi_ellipse = double(phi_ellipse);
    phi_ellipse( phi_ellipse == 1 ) = -2;
    phi_ellipse( phi_ellipse == 0 ) = 2;
    phi_ellipse = double((phi_ellipse > 0).*(bwdist(phi_ellipse < 0)- 0.5) - (phi_ellipse < 0).*(bwdist(phi_ellipse > 0)- 0.5));
    
    %=======================
    %       Edge term
    %=======================
    edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
    
    %=======================
    %      Update phi
    %=======================
    phi=phi + timestep*(mu*distRegTerm ...
        + lambda*edgeTerm ...           % Edge term
        + alfa*areaTerm ...             % Area term
        + gamma * phi_ellipse);         % Ellipse term

    phi = double((phi > 0).*(bwdist(phi < 0)- 0.5) - (phi < 0).*(bwdist(phi > 0)- 0.5));  
end


function f = phiArea(phi)
% count the number of pixels inside phi
phiBw = ~im2bw(phi);
phiStats = regionprops(phiBw, 'Area', 'PixelIdxList');
phiLabel = bwlabel(phiBw, 8);
phiIdxInLabel = find(phiLabel == 1);
for i = 1:length(phiStats)
    isPhi = isequal(phiStats(i,1).PixelIdxList(:), phiIdxInLabel);
    if isPhi == 1
        f = phiStats(i,1).Area;
        break;
    end
end


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