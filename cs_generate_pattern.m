function M = cs_generate_pattern (dim, accel, q, seed, pf)
%%CS_GENERATE_PATTERN  Random binary phase encode mask.
if nargin < 3, q = 1; end
if nargin < 4, seed = 11235813; end
if nargin < 5, pf = 1; end
if is_octave()
  rand('state', seed);
else
  rng(seed);
end
if numel(dim) == 2   % 2D mask
  nro = dim(1);
  npe = dim(2);
  nacq = round(npe / accel);  % requested # of lines in acquisition
  P = mask_pdf_1d(npe, nacq, q, pf)';
  C =  [];  % hackity-hack
  while true
    M = rand(npe,1) <= P;
    if sum(M) == nacq, break; end
  end
  % remove partial Fourier plane and compensate sampling density
  M = logical(M);
  M = repmat(M, [1 nro]);
elseif numel(dim) == 3  % 3D mask
  nro = dim(2);
  npe1 = dim(1);
  npe2 = dim(3);
  nacq = round(npe1*npe2 / accel);
  [P C] = mask_pdf_2d([npe1 npe2], nacq, q, pf);
  R = rand(npe1, npe2);
  M = R <= P;
  nchosen = sum(M(:));
  if nchosen > nacq   % correct for inexact number chosen
    outerOn = find(M & P ~= 1);
    numToFlip = nchosen - nacq;
    idxs = randperm(numel(outerOn));
    M(outerOn(idxs(1:numToFlip))) = false;
  elseif nchosen < nacq
    outerOff = find(~M);
    idxs = randperm(numel(outerOff));
    numToFlip = nacq - nchosen;
    M(outerOff(idxs(1:numToFlip))) = true;
  end
  % ky x kz > kz x ky > kz x ky x kx > ky x kx x kz
  M = shiftdim(repmat(shiftdim(M, 1), [1 1 nro]), 1);
  M = ifftshift(M);
  if nargout > 1
    C = shiftdim(repmat(shiftdim(C, 1), [1 1 nro]), 1);
    C = ifftshift(C);
  end
else
  error('Pattern dimension must be 2 or 3.');
end


function P = mask_pdf_1d (n, norm, q, pf)
%%VUCSMASKPDF Symmetric array of sampling probabilities.
if nargin < 4, pf = 1; end
ks = (1:n) - ceil(n/2) - 1;
kmax = floor(n / 2);
npf = round(pf*n);
klo = ks(n - npf + 1);
for kw = 1:kmax
  P = pdf(ks, kw, klo, q);
  if sum(P) >= norm, break; end
end
P = fftshift(P);
if mod(n, 2), P = [1; P]; end


function p = pdf (k, kw, klo, q)
p = (abs(k) / kw).^(-q);
p(k == 0) = 0;
p(abs(k) <= kw) = 1;
p(k < klo) = 0;


function [P C] = mask_pdf_2d (dims, norm, q, sense_factor)
if nargin < 4
  sense_factor = 1;
  if nargin < 3, q = 1; end
end
nz = dims(2);
ny = dims(1);
yc = round(ny/2);
zc = round(nz/2);
rmax = sqrt((ny-yc)^2 + (nz-zc)^2);
[Z Y] = meshgrid(1:nz, 1:ny);
R = sqrt((Y-yc).^2 + (Z-zc).^2);
Z = abs(Z - nz/2 - 0.5);
Y = abs(Y - ny/2 - 0.5);
for rw = 1:round(rmax)
  P = ones(ny, nz) / sense_factor;
  C = Z <= rw & Y <= rw;
  W = Z > rw | Y > rw;
  P(W) = (R(W) / rw) .^ (-q);
  if sum(P(:)) >= norm, break; end
end
if sense_factor > 1
  if round(sense_factor) ~= sense_factor
    error('Only integer SENSE factors have been implemented.');
  end
  D = C;
  for s = 2:sense_factor
    D(s:sense_factor:end, :) = false;
  end
  P(C) = 0;
  P(D) = 1;
end


function r = is_octave ()
persistent x
if isempty(x)
  x = exist('OCTAVE_VERSION', 'builtin');
end
r = x;
  
