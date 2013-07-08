function img = cs_recon_cartesian (data, mu)
%%CS_RECON_CARTESIAN  Spatial-only 2D and 3D CS recon.
%
% Perform a total variation constrained compressed sensing reconstruction
% of undersampled Cartesian data.  Samples not acquired are expected to be
% zeros in the arrays.
t = tic;
ninner = 5;
nbreg = 10;
lambda = 4.0; 
if nargin < 2; mu = 20; end
gamma = 1;
if ndims(data) == 2
  img = cs_recon_cartesian_2d(data, lambda, nbreg, ninner, mu, gamma);
elseif ndims(data) == 3
  img = cs_recon_cartesian_3d(data, lambda, nbreg, ninner, mu, gamma);
else
  data_size = size(data);
  data = reshape(data,size(data,1),size(data,2),size(data,3),[]);
  img = zeros(size(data),'single');
  for k = 1:size(data,4)
    img(:,:,:,k) = cs_recon_cartesian_3d(data(:,:,:,k), lambda, ...
      nbreg, ninner, gamma);
  end
  img = reshape(img, data_size);
end
img = abs(img);
toc(t);


function u = cs_recon_cartesian_2d(f, lambda, nbreg, ninner, mu, gamma)
[rows,cols] = size(f);
R = f ~= 0.0;   % not perfect, but good enough
% normalize the data so that standard parameter values work
norm_factor = get_norm_factor(R,f);
% perform recon in single precision to save RAM, since SNRs do not warrant
% a double precision recon anyway.
R = single(R);
f = single(f);
f = norm_factor * f;
% Reserve memory for the auxillary variables
f0 = f;
u = zeros(rows,cols,'single');
x = zeros(rows,cols,'single');
bx = zeros(rows,cols,'single');
y = zeros(rows,cols,'single');
by = zeros(rows,cols,'single');
% Build Kernels
scale = sqrt(rows*cols);
murf = ifft2(mu * R .* f) * scale;
uker = zeros(rows,cols,'single');
uker(1,1) = 4;
uker(1,2) = -1;
uker(2,1) = -1;
uker(rows,1) = -1;
uker(1,cols) = -1;
uker = 1 ./ (mu * R + lambda * fft2(uker) + gamma);
%  Do the reconstruction
for outer = 1:nbreg
  for inner = 1:ninner
    % update u
    %fprintf('i:%d, o:%d\n', inner, outer);
    rhs = murf + lambda*Dxt(x-bx) + lambda*Dyt(y-by) + gamma*u;
    u = ifft2(fft2(rhs) .* uker);
    % update x and y
    ax = Dx(u) + bx;
    ay = Dy(u) + by;
    [x y] = shrink2(ax, ay, 1/lambda);
    % update bregman parameters
    bx = ax - x;
    by = ay - y;
  end
  f = f + f0 - R .* fft2(u) / scale;
  murf = ifft2(mu * R .* f) * scale;
end
% undo the normalization so that results are scaled properly
u = u / norm_factor / scale;


function U = cs_recon_cartesian_3d(F, lambda, nbreg, ninner, mu, gamma)
%%VUCSRECONSB3D 3-D TV-regularized partial Fourier recon.
M = F ~= 0.0;
[rows cols slices] = size(F);
% normalize the data so that standard parameter values work
norm_factor = get_norm_factor(F,M);
F = single(F);
F = norm_factor * F;
% Reserve memory for the auxillary variables
F0 = F;
U = zeros(rows,cols,slices,'single');
X = zeros(rows,cols,slices,'single');
Bx = zeros(rows,cols,slices,'single');
Y = zeros(rows,cols,slices,'single');
By = zeros(rows,cols,slices,'single');
Z = zeros(rows,cols,slices,'single');
Bz = zeros(rows,cols,slices,'single');
% Build Kernels
scale = sqrt(rows*cols);   % TODO: should there by a slices factor here?
muMF = ifftn(mu * M .* F) * scale;
K = zeros(rows,cols,slices,'single');
K(1,1,1) = 8;
K(2,1,1) = -1;
K(1,2,1) = -1;
K(1,1,2) = -1;
K(rows,1,1) = -1;
K(1,cols,1) = -1;
K(1,1,slices) = -1;
K = 1 ./ (mu * M + lambda * fftn(K) + gamma);
nu = 1;
%  Do the reconstruction
for outer = 1:nbreg
  for inner = 1:ninner
    % update u
    RHS = muMF + lambda*Dxt(X-Bx) + lambda*Dyt(Y-By) + ...
      lambda*Dzt(Z-Bz)/nu + gamma*U;
    U = ifftn(fftn(RHS) .* K);
    % update x and y
    Ax = Dx(U) + Bx;
    Ay = Dy(U) + By;
    Az = nu * Dz(U) + Bz;
    [X Y Z] = shrink3(Ax, Ay, Az, 1/lambda);  
    % update bregman parameters
    Bx = Ax - X;
    By = Ay - Y;
    Bz = Az - Z; 
  end
  F = F + F0 - M .* fftn(U) / scale;
  muMF = ifftn(mu * M .* F) * scale;
end
% undo the normalization so that results are scaled properly
U = single(U / norm_factor / scale);


function norm_factor = get_norm_factor(R,f)
norm_factor = 1/norm(f(:)/size(R==1,1));


function d = Dx(u)
n = size(u,1);
d = u - u([n 1:n-1],:,:);

function d = Dxt(u)
d = u - u([2:size(u,1) 1],:,:);

function d = Dy(u)
n = size(u,2);
d = u - u(:,[n 1:n-1],:);

function d = Dyt(u)
d = u - u(:,[2:size(u,2) 1],:);  

function d = Dz(u)
n = size(u,3);
d = u - u(:,:,[n 1:n-1]);

function d = Dzt(u)
d = u - u(:,:,[2:size(u,3) 1]); 

function xs = shrink(x, gamma)
p = 1;
s = abs(x);
t = gamma ./ (s .^ (1 - p));
%t = gamma ./ sqrt(s);
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
xs = ss .* x;

function [xs, ys] = shrink2(x, y, lambda)
p = 1;
s = sqrt(x.*conj(x) + y.*conj(y));
t = lambda ./ (s .^ (1 - p));
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
xs = ss .* x;
ys = ss .* y;

function [xs, ys, zs] = shrink3(x, y, z, lambda)
p = 1;
s = sqrt(x.*conj(x) + y.*conj(y) + z.*conj(z));
t = lambda ./ (s .^ (1 - p));
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
xs = ss .* x;
ys = ss .* y;
zs = ss .* z;

