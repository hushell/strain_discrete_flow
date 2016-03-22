addpath(genpath('external/'));
addpath('tools');

img_path = 'data/';
img_name = 'initial_cropped.png';
res_path = 'results/navier_stokes/';

% read a reference image
im1=imread([img_path img_name]);
im1=im2double(im1(1:200,1:200));
[m,n,~] = size(im1);

saveFile = [res_path '/' img_name '_NavierStokes_GT.mat'];
if ~exist(saveFile, 'file')
    % initial random displacement field
    options.bound = 'per';
    v = perform_blurring(randn(m,n,2), 40, options);
    [tmp,v] = compute_hodge_decompositon(v,options);
    v = perform_vf_normalization(v);

    % Navier-Stokes dynamic system
    fprintf('Running Navier-Stokes dynamic system.\n');
    dt = .3;
    options.viscosity = 100*dt; % diffusion per frame
    options.advspeed = 1*dt; % advection per frame
    options.viscosity_texture = .3*dt; % diffusion of the texture
    options.texture_histo = 'linear'; % fix the contrast
    options.display = 0;
    options.niter_fluid = 100;

    % solve the PDE: vlist = {v for all iterations}, A = {all the warpped images}
    [vlist,A] = perform_fluid_dynamics(v,im1,options);

    % display some of A
    %sel = round( linspace(1,size(A,3),6) );
    %B = mat2cell(A(:,:,sel),n,n,ones(6,1));
    %clf;
    %imageplot(B);
    
    save(saveFile, 'vlist', 'A');
else
    load(saveFile);
end

v20 = vlist{20}*20;
vx = v20(:,:,1);
vy = v20(:,:,2);
len_v = max(abs(vx(:)), abs(vy(:)));
mean_v = mean(len_v);

vyWarp = vy(end:-1:1,:);
vxWarp = vx(end:-1:1,:);
im2=warpImage(im1,vxWarp,vyWarp); % transform im1 by the vector field [vx,vy]

%% Optical flow 
of_params.alpha=2*255;
of_params.d=40*255;
of_params.eta=0.005*255;
of_params.nlevels=4;
of_params.wsize=2;
%of_params.topwsize=10;
of_params.topwsize=2*mean_v;
of_params.nTopIterations = 60;
of_params.nIterations= 30;
of_params.useSIFT = 1;
of_params.cellSize = 3;
of_params.gridSpacing = 1;
of_params.patchSize = 12;

[ux,uy,energylist] = opt_flow_pair(im2, im1, of_params);

figure; cquiver(vx(1:5:end,1:5:end),vy(1:5:end,1:5:end))
figure; cquiver(ux(1:5:end,1:5:end),uy(1:5:end,1:5:end))

% Error analysis
len_v = sqrt((vx(:)).^2 + (vy(:)).^2);
len_u = sqrt((ux(:)).^2 + (uy(:)).^2);
anal.mean_v = mean(len_v);
anal.std_v = std(len_v);

diff_len = abs(1 - len_u ./ (len_v+1e-20));
anal.error_len = mean(diff_len(:));
%anal.std_len = std(abs((len_u-len_v))./(len_v+1e-20));
anal.std_len = std(diff_len);

diff_ang = vx(:).*ux(:) + vy(:).*uy(:);
diff_ang = acos( diff_ang ./ ((len_u+1e-20) .* (len_v+1e-20)) );
anal.error_ang = mean(diff_ang) / pi * 180;
