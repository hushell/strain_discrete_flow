function [u, v, w] = tv_l1_optical_flow(I1, I2, options) 
% Inputs: I1, I2 are two consecutive images (only grayscale allowed)
% Outputs: (u,v) is the flow, where size(u) = size(v) = size(I1)
%          w is the estimated illumination changes
% Author: Thomas Pock and Chen Yunjin and Shell Hu

lambda = options.lambda;
beta = options.beta;
warps = options.warps;
max_iter = options.max_iter;
pyramid_levels = options.pyramid_levels;
pyramid_factor = options.pyramid_factor;
check = options.check;

[M N C] = size(I1);

num_dual_vars = 6;

width_Pyramid = cell(pyramid_levels,1);
height_Pyramid = cell(pyramid_levels,1);

I1_Pyramid = cell(pyramid_levels,C);
I2_Pyramid = cell(pyramid_levels,C);

% precompute image sizes
width_Pyramid{1} = N;
height_Pyramid{1} = M;
for i = 2:pyramid_levels
  width_Pyramid{i} = pyramid_factor*width_Pyramid{i-1};
  height_Pyramid{i} = pyramid_factor*height_Pyramid{i-1};
  if min(width_Pyramid{i}, height_Pyramid{i}) < 16
    pyramid_levels = i;
    break;
  end
end
for i = 1:pyramid_levels
  width_Pyramid{i} = round(width_Pyramid{i});
  height_Pyramid{i} = round(height_Pyramid{i});
end

% set up image pyramides
for i = 1:pyramid_levels
  if i == 1
    for j=1:C
      I1_Pyramid{1,j} = I1(:,:,j);
      I2_Pyramid{1,j} = I2(:,:,j);
    end
  else
    for j=1:C    
      I1_Pyramid{i,j} = imresize(I1_Pyramid{i-1,j}, ...
                                 [height_Pyramid{i} width_Pyramid{i}], 'bicubic');
      I2_Pyramid{i,j} = imresize(I2_Pyramid{i-1,j}, ...
                                 [height_Pyramid{i} width_Pyramid{i}], 'bicubic');     
    end
  end
end

% main loop
for level = pyramid_levels:-1:1;
  scale = pyramid_factor^(level-1);
  
  M = height_Pyramid{level};
  N = width_Pyramid{level};
 
  if level == pyramid_levels
    % initialization  
    u = zeros(M,N);
    v = zeros(M,N);
    w = zeros(M,N);
    p = zeros(M,N,num_dual_vars);   
  else
    rescale_factor_u = width_Pyramid{level+1}/width_Pyramid{level};
    rescale_factor_v = height_Pyramid{level+1}/height_Pyramid{level};
    
    % prolongate to finer grid    
    u = imresize(u,[M N], 'nearest')/rescale_factor_u;    
    v = imresize(v,[M N], 'nearest')/rescale_factor_v;
    w = imresize(w, [M N], 'nearest');
    
    p_tmp = p;
    p = zeros(M,N,num_dual_vars); 
    for i=1:num_dual_vars
      p(:,:,i) = imresize(p_tmp(:,:,i),[M N], 'nearest');
    end
  end
  
  I1 = zeros(M,N,C);
  I2 = zeros(M,N,C);
  for j=1:C
    I1(:,:,j) = I1_Pyramid{level,j};
    I2(:,:,j) = I2_Pyramid{level,j};
  end
  
  fprintf('*** level = %d\n', level);
  
  [u, v, w, p] = champolle_pock_primal_dual_alg(I1, I2, ...
                        u, v, w, p, lambda, warps, max_iter, scale, beta, check);
  
  fprintf('=== Processing...%.2f/100\n', (pyramid_levels + 1 - level)/pyramid_levels*100);
  drawnow;
end
fprintf('Done!\n');
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, v, w, p] = champolle_pock_primal_dual_alg(I1, I2, ...
                               u, v, w, p, lambda, warps, maxits, scale, beta, check)
% [M N C] = size(I1);
[n_row,n_col] = size(I1);
N             = n_row * n_col;
nabla         = make_nabla(n_col,n_row); % Kxf, Kyf for gradient
num_dual_vars = 3;

% stepwidth
L = sqrt(8);

tau = 1/L;
sigma = 1/L;

% some parameters
epsilon_u = 0;
epsilon_w = 0;
% vectorization
vector_u   = reshape(u',N,1); % primal
vector_v   = reshape(v',N,1); % primal
vector_w   = reshape(w',N,1); % auxilary primal illumination
vector_u_  = vector_u; % primal_bar
vector_v_  = vector_v; % primal_bar
vector_w_  = vector_w; % primal_bar
vector_p   = zeros(2*N,num_dual_vars);
for i = 1:num_dual_vars;
  vector_p(1:N,i)     = reshape(p(:,:,2*i-1)',N,1);
  vector_p(N+1:2*N,i) = reshape(p(:,:,2*i)',N,1);
end
nabla_t = nabla'; % K^T for divergence

% inner loops for warping
for j = 1:warps
  vector_u0 = vector_u;
  vector_v0 = vector_v;
  u0 = (reshape(vector_u0,n_col,n_row))';
  v0 = (reshape(vector_v0,n_col,n_row))';
  
  % warping
  fprintf('tv-l1-of-pd: warp = %d\n', j);
  
  [I_x, I_y, I_t, I2_warped] = warping(I1, I2, vector_u0, vector_v0);  
  % vectorization
  vector_I_x = reshape(I_x',N,1);
  vector_I_y = reshape(I_y',N,1);
  vector_I_t = reshape(I_t',N,1);
  I_grad_sqr = max(1e-09, vector_I_x.^2 + vector_I_y.^2 + beta*beta);
  
  % iterations of Champolle-Pock alg
  for k = 0:maxits-1
    %% DUAL
    % compute derivatives and update dual variable: p_tilde = p + sigma*K*uv_bar
    vector_p(:,1:2) = (vector_p(:,1:2) + sigma*nabla*[vector_u_,vector_v_])/(1 + sigma*epsilon_u);
    vector_p(:,3)   = (vector_p(:,3)   + sigma*nabla*vector_w_)/(1 + sigma*epsilon_w); 
    
    % reprojection to |pu| <= 1
    norm = repmat(sqrt(vector_p(1:N,1).^2 + vector_p(N+1:2*N,1).^2 + ...
        vector_p(1:N,2).^2 + vector_p(N+1:2*N,2).^2),[2,2]);
    reprojection = max(1.0,norm);
    vector_p(:,1:2) = vector_p(:,1:2)./reprojection;
    
    norm = repmat(sqrt(vector_p(1:N,3).^2 + vector_p(N+1:2*N,3).^2),[2,1]);
    reprojection = max(1.0,norm);
    vector_p(:,3) = vector_p(:,3)./reprojection;
    
    %% PRIMAL
    % remember old u,v,w
    vector_u_ = vector_u;
    vector_v_ = vector_v;
    vector_w_ = vector_w;
  
    % uv_tilde = uv - tau*K^T*p
    vector_u = vector_u - tau * nabla_t * vector_p(:,1);
    vector_v = vector_v - tau * nabla_t * vector_p(:,2);
    vector_w = vector_w - tau * nabla_t * vector_p(:,3);
    
    % prox operator for u,v,w
    rho = vector_I_t + (vector_u - vector_u0).*vector_I_x + ...
        (vector_v - vector_v0).*vector_I_y + beta*vector_w;
    idx1 = rho      < - tau*lambda*I_grad_sqr;
    idx2 = rho      >   tau*lambda*I_grad_sqr;
    idx3 = abs(rho) <=  tau*lambda*I_grad_sqr;
    
    vector_u(idx1) = vector_u(idx1) + tau*lambda*vector_I_x(idx1);
    vector_v(idx1) = vector_v(idx1) + tau*lambda*vector_I_y(idx1);
    vector_w(idx1) = vector_w(idx1) + tau*lambda*beta;
    
    vector_u(idx2) = vector_u(idx2) - tau*lambda*vector_I_x(idx2);
    vector_v(idx2) = vector_v(idx2) - tau*lambda*vector_I_y(idx2);
    vector_w(idx2) = vector_w(idx2) - tau*lambda*beta;
    
    vector_u(idx3) = vector_u(idx3) - rho(idx3).*vector_I_x(idx3)./I_grad_sqr(idx3);
    vector_v(idx3) = vector_v(idx3) - rho(idx3).*vector_I_y(idx3)./I_grad_sqr(idx3);
    vector_w(idx3) = vector_w(idx3) - rho(idx3).*beta./I_grad_sqr(idx3);
    
    vector_u_ = 2*vector_u - vector_u_;
    vector_v_ = 2*vector_v - vector_v_;
    vector_w_ = 2*vector_w - vector_w_;
    
    if mod(k,check) == 0
      u = (reshape(vector_u,n_col,n_row))';
      v = (reshape(vector_v,n_col,n_row))';
      w = (reshape(vector_w,n_col,n_row))';      
      figure(109); show_flow(u,v,beta*w,I1,I2_warped + (u-u0).*I_x + (v-v0).*I_y + beta*w);
      fprintf('tv-l1-motion-primal-dual: it = %d\n', k)
    end

  end
  
  % recover matrix representation
  u = (reshape(vector_u,n_col,n_row))';
  v = (reshape(vector_v,n_col,n_row))';
  w = (reshape(vector_w,n_col,n_row))';
  for i = 1:num_dual_vars;
      p(:,:,2*i-1) = (reshape(vector_p(1:N,i),n_col,n_row))';
      p(:,:,2*i)   = (reshape(vector_p(N+1:2*N,i),n_col,n_row))';
  end
  % filter strong outliers
  u = peakfilt(u);
  v = peakfilt(v);
end

%----------------------------------------
function [nabla] = make_nabla(M,N)

% Kxf{matrix for gradient_x}
row  = zeros(1,M*N*2);
col = zeros(1,M*N*2);
val  = zeros(1,M*N*2);

cnt = 1;

for y=1:M-1
  for x=1:N
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x-1)*M;
    val(cnt) = -1;
    cnt = cnt+1;
    
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x-1)*M+1;
    val(cnt) = 1;
    cnt = cnt+1;
  end
end
row = row(1:cnt-1);
col = col(1:cnt-1);
val = val(1:cnt-1);

Kxf = sparse(row,col,val,M*N,M*N);

% Kyf(matrix for gradient_y)
row  = zeros(1,M*N*2);
col = zeros(1,M*N*2);
val  = zeros(1,M*N*2);

cnt = 1;

for y=1:M
  for x=1:N-1
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x-1)*M;
    val(cnt) = -1;
    cnt = cnt+1;
    
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x)*M;
    val(cnt) = 1;
    cnt = cnt+1;
  end
end
row = row(1:cnt-1);
col = col(1:cnt-1);
val = val(1:cnt-1);

Kyf = sparse(row,col,val,M*N,M*N);

nabla = [Kxf;Kyf];

%----------------------------------------
function [I_x, I_y, I_t, I2_warped] = warping(I1, I2, u, v)

[M N C] = size(I1);
u = (reshape(u,N,M))';
v = (reshape(v,N,M))';

idx = repmat([1:N], M,1);
idy = repmat([1:M]',1,N);
 
idxx = idx + u;
idyy = idy + v;
m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);

idxx = max(1,min(N,idxx));
idxm = max(1,min(N,idxx-0.5));
idxp = max(1,min(N,idxx+0.5));

idyy = max(1,min(M,idyy));
idym = max(1,min(M,idyy-0.5));
idyp = max(1,min(M,idyy+0.5));

I2_warped = interp2(I2,idxx,idyy,'cubic');
I2_x_warped = interp2(I2,idxp,idyy,'cubic') - interp2(I2,idxm,idyy,'cubic');
I2_y_warped = interp2(I2,idxx,idyp,'cubic') - interp2(I2,idxx,idym,'cubic');

% use everage to improve accuracy
I_x = I2_x_warped;
I_y = I2_y_warped;

I_t = I2_warped  - I1;

% boundary handling
I_x(m) = 0.0;
I_y(m) = 0.0;
I_t(m) = 0.0;

%----------------------------------------
function [u] = peakfilt(u);
u_ = medfilt2(u,[3 3],'symmetric');
diff = abs(u-u_);
v = mean(abs(u_(:)));
%sum(sum(diff > v))/size(u,1)/size(u,2)
u(diff > v) = u_(diff > v);

%----------------------------------------
function show_flow(u,v,c,I1,I2)
[M N] = size(u);
% find robust max flow for better visualization
magnitude = (u.^2 + v.^2).^0.5;  
max_flow = prctile(magnitude(:),95);

tmp = zeros(M,N,2);
tmp(:,:,1) = min(max(u,-max_flow),max_flow);
tmp(:,:,2) = min(max(v,-max_flow),max_flow);
if max(tmp(:)) ~= min(tmp(:))  
  subplot(2,2,1); imshow(I1,[0 1]); title('I1');
  subplot(2,2,2); imshow(I2,[0 1]); title('warped I2 -> I1');
  subplot(2,2,3); imshow(uint8(flowToColor(tmp)),[]); title('flow');
  %subplot(2,2,3), imshow(sqrt(tmp(:,:,1).^2 + tmp(:,:,2).^2),[]);  
  %subplot(2,2,4), imshow(c,[-0.01 0.01]); drawnow;
  subplot(2,2,4); imshow(c,[]); title('illumination changes');
  drawnow;  
end

