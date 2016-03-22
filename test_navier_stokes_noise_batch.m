function test_navier_stokes_noise_batch(imFile, savePath, paraName, valRange, useSIFT)

if nargin < 5
    useSIFT = 1;
end

[~,imName] = fileparts(imFile);

% read a reference image
im1=imread(imFile);
im1=im2double(im1);
[m,n,~] = size(im1);

saveFile = [savePath '/' imName '_NavierStokes_GT.mat'];
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

%% We choose the v from iter 20
fprintf('Running opt_flow experiments.\n');
for s = 5:5:50 % s: scalar controls the length of the vectors
    fprintf('scalar = %d\n', s);
    v20 = vlist{20}*s;
    vx = v20(:,:,1);
    vy = v20(:,:,2);
    len_v = max(abs(vx(:)), abs(vy(:)));
    mean_v = mean(len_v);

    vyWarp = vy(end:-1:1,:);
    vxWarp = vx(end:-1:1,:);
    im2=warpImage(im1,vxWarp,vyWarp); % transform im1 by the vector field [vx,vy]
    
    nv = var(im2(:));
    for t = [0 nv/2] % t: Gaussian noise var
        fprintf('var = %f\n', t);
        im2 = imnoise(im2, 'gaussian', 0, t);

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
        of_params.useSIFT = useSIFT;
        of_params.cellSize = 3;
        of_params.gridSpacing = 1;
        of_params.patchSize = 12;

        for i = 1:length(valRange)
            fprintf('%s = %f\n', paraName, valRange(i));
            of_params.(paraName) = valRange(i);
            [ux,uy,energylist] = opt_flow_pair(im2, im1, of_params);
            save([savePath '/' imName '_' paraName num2str(valRange(i)) '_s' num2str(s) '_var' num2str(t) '.mat'], ...
                'ux', 'uy', 'energylist');
        end
    end
end
