function [d] = lbfgsProd(g,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag)
% BFGS Search Direction
%
% This function returns the (L-BFGS) approximate inverse Hessian,
% multiplied by the negative gradient

% Determine if we are working on the GPU
if strcmp(class(S), 'gpuArray')
	gpu = 1;
else
    gpu = 0;
end

% Set up indexing
[nVars,maxCorrections] = size(S);
if lbfgs_start == 1
	ind = 1:lbfgs_end;
	nCor = lbfgs_end-lbfgs_start+1;
else
	ind = [lbfgs_start:maxCorrections 1:lbfgs_end];
	nCor = maxCorrections;
end

if gpu
    al = gpuArray.zeros(nCor,1);
    be = gpuArray.zeros(nCor,1);
else
    al = zeros(nCor,1);
    be = zeros(nCor,1);
end

d = -g;
for j = 1:length(ind)
	i = ind(end-j+1);
	al(i) = (S(:,i)'*d)/YS(i);
	d = d-al(i)*Y(:,i);
end

% Multiply by Initial Hessian
d = Hdiag*d;

for i = ind
	be(i) = (Y(:,i)'*d)/YS(i);
	d = d + S(:,i)*(al(i)-be(i));
end
