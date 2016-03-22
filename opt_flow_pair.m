function [vx,vy,energylist] = opt_flow_pair(imfile1, imfile2, of_params)

if isstr(imfile1) 
    assert(isstr(imfile2));
    im1=imread(imfile1);
    im2=imread(imfile2);
else
    im1 = imfile1;
    im2 = imfile2;
end

if nargin < 3
    of_params.alpha=2*255; % smooth coeff
    of_params.d=40*255; % smooth thres
    of_params.eta=0.005*255; % small displacement coeff
    of_params.nlevels=4; % n levels of pyramid
    of_params.wsize=2; % search window size below top
    of_params.topwsize=10; % search window size of top level
    of_params.nTopIterations = 60; % n iters of top level
    of_params.nIterations= 30; % n iters below top
    of_params.useSIFT = 1; % whether use SIFT features
    of_params.cellSize = 3; % SIFT cell size
    of_params.gridSpacing = 1; % SIFT grid spacing
    of_params.patchSize = 12; % SIFT patch size = 4*cellSize x 4*cellSize
end

if of_params.useSIFT
    im1=im2double(im1);
    im2=im2double(im2);
    cellsize = of_params.cellSize;
    gridspacing = of_params.gridSpacing;
    feat1 = mexDenseSIFT(im1,cellsize,gridspacing);
    feat2 = mexDenseSIFT(im2,cellsize,gridspacing);
else
    im1 = uint8(im1);
    im2 = uint8(im2);
    patchsize = of_params.patchSize;
    feat1 = topLeftPatchExtract(im1,patchsize);
    feat2 = topLeftPatchExtract(im2,patchsize);
end

tic;
[vx,vy,energylist] = opt_flow(feat1,feat2,of_params);
toc
vy = vy(end:-1:1,:);
vx = vx(end:-1:1,:);

