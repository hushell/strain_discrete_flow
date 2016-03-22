function results = eval_navier_stokes_noise_batch(imFile, savePath, paraName, valRange)

[~,imName] = fileparts(imFile);

% read a reference image
im1=imread(imFile);
im1=im2double(im1);
[m,n,~] = size(im1);

saveFile = [savePath '/' imName '_NavierStokes_GT.mat'];
load(saveFile);

results = cell(length(5:5:50), 2, length(valRange));
fprintf('Evaluating opt_flow experiments.\n');
ss = 0;
for s = 5:5:50
    ss = ss + 1;
    fprintf('s = %d\n', s);
    %s = 20; % scalar controls the length of the vectors
    v20 = vlist{20}*s;
    vx = v20(:,:,1);
    vy = v20(:,:,2);

    % im2 = warp im1 by v20
    vyWarp = vy(end:-1:1,:);
    vxWarp = vx(end:-1:1,:);
    im2=warpImage(im1,vxWarp,vyWarp);
    
    nv = var(im2(:));
    
    tt = 0;
    for t = [0 nv/2]
        tt = tt + 1;
        
        for i = 1:length(valRange)
            fprintf('%s = %f\n', paraName, valRange(i));
            saveFile = [savePath '/' imName '_' paraName '_' num2str(valRange(i)) 's_' num2str(s) 'noisvar_' num2str(t) '.mat'];
            if exist(saveFile,'file')
                load(saveFile);
                anal = error_anal(vx,vy,ux,uy);
                results{ss,tt,i} = anal;
            end
        end
    end
end



function anal = error_anal(vx,vy,ux,uy)
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