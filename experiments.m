addpath(genpath('external/'));
addpath('tools');

img_path = 'data/';
img_name = 'initial_cropped.png';
res_path = 'results/navier_stokes/';

% try diff patch sizes
useSIFT = 0;
test_navier_stokes_noise_batch([img_path img_name], res_path, 'patchSize', [1:4]*4, useSIFT);
results = eval_navier_stokes_noise_batch([img_path img_name], res_path, 'patchSize', [1:4]*4);

%% Some analyses
% 10 scales x 2 noise levels
len_s = 10;
len_noise = 2;
min_len = zeros(len_s,len_noise);
min_ang = zeros(len_s,len_noise);
min_stdlen = zeros(len_s,len_noise);

wc_len = min_len;
wc_ang = min_ang;
wc_slen = min_len;
mean_v = zeros(1,len_s);
for i = 1:len_s
    for j = 1:len_noise
    fprintf('s = %d\n', i);
    si = results(i,j,:);
    [min_len(i,j),wc_len(i,j)] = min(cellfun(@(r) r.error_len, si));
    [min_ang(i,j),wc_ang(i,j)] = min(cellfun(@(r) r.error_ang, si));
    [min_stdlen(i,j),wc_slen(i,j)] = min(cellfun(@(r) r.std_len, si));
    
    mean_v(i) = si{1,1}.mean_v;
    end
end

figure; plot(mean_v,wc_len(:,1),'DisplayName','wc_len','YDataSource','wc_len');
figure; plot(mean_v(:),min_len(:,1),mean_v(:),min_len(:,2),'DisplayName','min_len','YDataSource','min_len');
figure; plot(mean_v,wc_ang(:,1),'b',mean_v,wc_ang(:,2),'r','DisplayName','wc_ang','YDataSource','wc_ang');
figure; plot(mean_v,min_ang(:,1),mean_v,min_ang(:,2),'DisplayName','min_len','YDataSource','min_len');

load([img_path img_name '_patchSize4_s10_var0.mat');
%figure; cquiver(vx(1:5:end,1:5:end),vy(1:5:end,1:5:end))
figure; cquiver(ux(1:5:end,1:5:end),uy(1:5:end,1:5:end));

