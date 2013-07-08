clear all;

% EXAMPLE: 3D retrospective undersampling and recon of brain

% load image
load brains I256;
I = I256;

% generate a random sampling pattern
mask = cs_generate_pattern(size(I), 4);

% decimate the full Fourier data, leaving zeros in 'unacquired' spots
data = fft2(I) .* mask;

% run the recon
tic;
J = cs_recon_cartesian(single(data));
toc;

% view the results
figure(1);
subplot(121);
imagesc(I);
title('original');
colormap(gray);
subplot(122);
imagesc(J);
title('CS reconstruction');



%%
nz = 4;
data_3d = repmat(data, [1 1 nz]);
J_3d = cs_recon_cartesian(data_3d);

figure(2);
imagesc([I J_3d(:,:,1)]);
colormap(gray);
