function [ L,S ] = lowrank(weight2,oldx, i1, j1, t, f, l, s) 
[line, row] = size(weight2);
th = 0;
count = 1;
for I = 3:line-3
    for J = 3:row-3
        if (weight2(I,J)~=0 && weight2(I,J)== max(weight2(I-2:I+2,J)))
            Win = oldx(i1-t-1+I-f:i1-t-1+I+f,j1-t-1+J-f:j1-t-1+J+f);
            W = Win(:);
            Matrix(:,count) = W;
            sortW(count) = weight2(I,J); 
            count = count + 1;
        end
    end
end
clear W Win;
[Wsort, Index] = sort(sortW, 'DESCEND');
clear sortW sortW;
Matrix_new = zeros(size(Matrix));
if (length(Index) > t)
    Ind = Index(1:t);   
    Matrix_new = double(Matrix(:,Ind));
    clear Ind;
else
    Matrix_new = double(Matrix(:,Index));
    clear Index;
end
clear Matrix;
params.progTol = 1e-10;
params.optTol  = 1e-10;
params.MaxIter = 100;
params.store    = 1;
% Set parameters for split_SPCP
params.k   = 15; % rank bound on L 10
params.gpu = 0;  % set to 1 to run on GPU
params.lambdaL = l;%115  15  17
params.lambdaS = s;%0.8
[L,S] = solver_split_SPCP(Matrix_new,params);
end

