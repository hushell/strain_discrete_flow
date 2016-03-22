n = 256;
M = load_image('lena',n);
options.bound = 'per';
v = perform_blurring(randn(n,n,2), 40, options);
[tmp,v] = compute_hodge_decompositon(v,options);
v = perform_vf_normalization(v);

% iterative advection
dt = .7; niter = 30;
Ma = M; Msvg = {};
for i=1:niter
    Ma = perform_image_advection(Ma, dt*v, options);
    if mod(i,5)==1
        Msvg{end+1} = Ma;
    end
end
% display
clf;
imageplot(Msvg);

% Navier-Stokes
dt = .3;
options.viscosity = 100*dt; % diffusion per frame
options.advspeed = 1*dt; % advection per frame
options.viscosity_texture = .3*dt; % diffusion of the texture
options.texture_histo = 'linear'; % fix the contrast
options.display = 0;
options.niter_fluid = 100;
% solve the PDE
[vlist,A] = perform_fluid_dynamics(v,M,options);
% display
sel = round( linspace(1,size(A,3),6) );
B = mat2cell(A(:,:,sel),n,n,ones(6,1));
clf;
imageplot(B);

u = vlist{2};
figure; cquiver(u(1:5:end,1:5:end,1),u(1:5:end,1:5:end,2))