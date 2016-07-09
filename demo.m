% load images
%fpath = './frame10.png';
%fpath = './Mars-1.jpg';
fpath = './0275.png';
img_src1 = im2double(imread(fpath));
if size(img_src1,3) == 3
    img_src1 = rgb2gray(img_src1);
end
%fpath = './frame11.png';
%fpath = './Mars-2.jpg';
fpath = './0279.png';
img_src2 = im2double(imread(fpath));
if size(img_src2,3) == 3
    img_src2 = rgb2gray(img_src2);
end
%figure; imshow(img_src1);
%figure; imshow(img_src2);

% tv-l1 flow (coarse-to-fine)
options.lambda = 40;
options.beta   = 0.01;
max_iter = 50;
options.max_iter = round(max_iter);
check = 10;
options.check = round(check);
pyramid_levels = 1000;
options.pyramid_levels = round(pyramid_levels);
options.pyramid_factor = 0.9;
warps = 1;
options.warps = round(warps);
[u, v, illumination] = tv_l1_optical_flow(img_src1, img_src2, options);

%% post-processing
%% find robust max flow for better visualization
%magnitude = (u.^2 + v.^2).^0.5;  
%max_flow = prctile(magnitude(:),95);
%
%illumination = illumination - min(illumination(:));
%illumination = illumination/max(illumination(:));
%
%u = min(max(u,-max_flow),max_flow);
%v = min(max(v,-max_flow),max_flow);

%% write flow
%motfile = 'motion.png';
%illfile = 'illumination.png';
%[M N C] = size(img_src1);
%flow = zeros(M,N,2);
%flow(:,:,1) = u;
%flow(:,:,2) = v;
%imwrite(uint8(flowToColor(flow)),motfile);
%imwrite(illumination,illfile);

% append the estimated optical flow to the first frame to show the motion
% of every pixel
figure; imshow(img_src1);
hold on;
[n_row,n_col] = size(img_src1);
[X,Y] = meshgrid(1:10:n_col,1:10:n_row);
quiver(X,Y,u(1:10:end,1:10:end),v(1:10:end,1:10:end));
title('estimated optical flow appended to the first frame');

% warping
%sz = size(img_src1);
%[X,Y] = meshgrid(1:sz(2), 1:sz(1));
%back = interp2(img_src2,X+flow(:,:,1),Y+flow(:,:,2));
%figure; imagesc(abs(back-img_src1))
im1_back = warpImage(img_src2, u, v);
figure; imshow(im1_back); title('warped I2 -> I1');
