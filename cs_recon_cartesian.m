function img = cs_recon_cartesian (data, mu)
%%CS_RECON_CARTESIAN  Spatial-only 2D and 3D CS recon.
%
% Total variation (L1-norm of gradient) constrained compressed
% sensing reconstruction of undersampled Cartesian data.
% Samples not acquired are expected to be zeros in the arrays.
%
ninner = 5;
nbreg = 10;
lambda = 4.0;
if nargin < 2; mu = 20; end
gam = 1;
if ndims(data) == 2 || ndims(data) == 3
  img = cs_recon_cartesian_helper(data, lambda, nbreg, ninner, mu, gam);
else
  data_size = size(data);
  data = reshape(data, size(data,1), size(data,2), size(data,3), []);
  img = zeros(size(data),'single');
  for k = 1:size(data,4)
    img(:,:,:,k) = cs_recon_cartesian_helper(data(:,:,:,k), lambda, ...
      nbreg, ninner, gam);
  end
  img = reshape(img, data_size);
end
img = abs(img);


function img = cs_recon_cartesian_helper(data, lambda, nbreg, ninner, mu, gam)
data_size = size(data);
data_type = class(data);
data_ndims = ndims(data);
assert(data_ndims == 2 || data_ndims == 3);
mask = data ~= 0.0;   % not perfect, but good enough
% normalize the data so that standard parameter values work
norm_factor = get_norm_factor(mask, data);
data = norm_factor * data;
% Reserve memory for the auxillary variables
data0 = data;
img = zeros(data_size, data_type);
X = zeros([data_size, data_ndims], data_type);
B = zeros([data_size, data_ndims], data_type);
% Build Kernels
scale = sqrt(prod(data_size));
murf = ifftn(mu * mask .* data) * scale;
uker = zeros(data_size, data_type);
if data_ndims == 2
  uker(1,1) = 4;
  uker(1,2) = -1; uker(2,1) = -1;
  uker(end,1) = -1; uker(1,end) = -1;
else data_ndims == 3
  uker(1,1,1) = 8;
  uker(2,1,1) = -1; uker(1,2,1) = -1; uker(1,1,2) = -1;
  uker(end,1,1) = -1; uker(1,end,1) = -1; uker(1,1,end) = -1;
end
uker = 1 ./ (mu * mask + lambda * fftn(uker) + gam);
%  Do the reconstruction
for outer = 1:nbreg
  for inner = 1:ninner
    % update u
    %fprintf('i:%d, o:%d\n', inner, outer);
    rhs = murf + lambda*Dxyzt(X-B) + gam*img;
    img = ifftn(fftn(rhs) .* uker);
    % update x and y
    A = Dxyz(img) + B;
    X = shrink(A, 1/lambda);
    % update bregman parameters
    B = A - X;
  end
  data = data + data0 - mask .* fftn(img) / scale;
  murf = ifftn(mu * mask .* data) * scale;
end
% undo the normalization so that results are scaled properly
img = img / norm_factor / scale;


function norm_factor = get_norm_factor(mask, data)
norm_factor = 1/norm(data(:)/size(mask==1,1));

function D = Dxyz(u)
dx = Dx(u);
dy = Dy(u);
if ndims(u) == 3
  dz = Dz(u);
else
  dz = [];
end
D = cat(ndims(u)+1,dx,dy,dz);

function u = Dxyzt(X)
if ndims(X) == 4
  u = Dxt(X(:,:,:,1)) + Dyt(X(:,:,:,2)) + Dzt(X(:,:,:,3));
elseif ndims(X) == 3
  u = Dxt(X(:,:,1)) + Dyt(X(:,:,2));
end

function d = Dx(u)
n = size(u,1);
d = u - u([n 1:n-1], :, :);

function d = Dxt(u)
d = u - u([2:size(u,1) 1], :, :);

function d = Dy(u)
n = size(u,2);
d = u - u(:, [n 1:n-1], :);

function d = Dyt(u)
d = u - u(:, [2:size(u,2) 1], :);

function d = Dz(u)
n = size(u,3);
d = u - u(:, :, [n 1:n-1]);

function d = Dzt(u)
d = u - u(:, :, [2:size(u,3) 1]);

function Xs = shrink(X, gam)
p = 1;
s = abs(X);
t = gam ./ (s .^ (1 - p));
%t = gam ./ sqrt(s);
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
Xs = ss .* X;
