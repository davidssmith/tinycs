load brains I256;
I = I256;
%mask = vuCSRandMask(size(I),3);
mask = cs_generate_pattern(size(I), 3);
data = fft2(I) .* mask;
J = cs_recon_cartesian(data);

figure(1);
imagesc([I J]);
colormap(gray);


%%
nz = 8;
data_3d = repmat(data, [1 1 nz]);
J_3d = cs_recon_cartesian(data_3d);

figure(2);
imagesc([I J_3d(:,:,1)]);
colormap(gray);