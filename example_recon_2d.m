clear all; clf; clc;

% EXAMPLE: 2D retrospective undersampling and recon of brain

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

