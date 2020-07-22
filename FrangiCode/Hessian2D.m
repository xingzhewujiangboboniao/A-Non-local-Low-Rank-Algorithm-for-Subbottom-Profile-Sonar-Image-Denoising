function [Dxx,Dxy,Dyy] = Hessian2D(I,Sigma)
%% original Frangi Filter
if nargin < 2, Sigma = 1; end
[X,Y] = ndgrid(-round(3*Sigma):round(3*Sigma));
%%
%%%一阶导数
% DGaussx = (-1/(2*pi*Sigma^4)) * (X.*exp(-(X.^2 + Y.^2)/(2*Sigma^2)));
% DGaussy = (-1/(2*pi*Sigma^4)) * (Y.*exp(-(X.^2 + Y.^2)/(2*Sigma^2)));
%%

DGaussxx = 1/(2*pi*Sigma^4) * (X.^2/Sigma^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussxy = 1/(2*pi*Sigma^6) * (X .* Y)           .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussyy = DGaussxx';
Dxx = imfilter(I,DGaussxx,'conv');
Dxy = imfilter(I,DGaussxy,'conv');
Dyy = imfilter(I,DGaussyy,'conv');

%%
%%%先滤波再求导
% h = fspecial('gaussian',[11,11], 2);
% Id = imfilter(I, h);
% [line, row] = size(Id);
% 
% Dy=gradient(Id,'y');
% Dyy=(gradient(Dy,'y'));
% clear Dy;
% 
% Dx=gradient(Id,'x');
% Dxx=(gradient(Dx,'x'));
% Dxy=(gradient(Dx,'y'));
% clear Dx;
% % Y = [1  1  1  1  1  -2 -2 -2 -2 -2 1 1  1  1  1  ;
% %      1  1  1  1  1  -2 -2 -2 -2 -2 1 1  1  1  1  ;
% %      1  1  1  1  1  -2 -2 -2 -2 -2 1 1  1  1  1  ;
% %        ];
% % X = [1 1 1;
% %      1 1 1;
% %      1 1 1;
% %      1 1 1;
% %      1 1 1;
% %      -2 -2 -2;
% %      -2 -2 -2;
% %      -2 -2 -2;
% %      -2 -2 -2;
% %      -2 -2 -2;
% %      1   1  1;
% %      1   1  1;
% %      1   1  1;
% %      1   1  1;
% %       1  1  1;
% %            ];
% % Dx = imfilter(Id,X,'conv');
% % Dy = imfilter(Id,Y,'conv');
% % Dxx = imfilter(Dx,X,'conv');
% % Dyy = imfilter(Dy,Y,'conv');
% % Dxy = imfilter(Dx,Y,'conv');
% 
% Dx = zeros(line,row);
% Dy = zeros(line,row);
% Dxx = zeros(line,row);
% Dxy = zeros(line,row);
% Dyy = zeros(line,row);
% for i = 2:line-2
%     for j = 2:row-2
%         Dx(i,j) = (Id(i+1,j) - Id(i-1,j))/2;
%         Dy(i,j) = (Id(i,j+1) - Id(i,j-1))/2; 
%     end
% end
% for i = 2:line-1
%     for j = 2:row-1
%         Dxx(i,j) = (Dx(i+1,j) - Dx(i-1,j))/2;
%         Dyy(i,j) = (Dy(i,j+1) - Dy(i,j-1))/2; 
%         Dxy(i,j) = (Dx(i,j+1) - Dx(i,j-1))/2;
%     end
% end
%%
%%%显示二阶导数
% mag = sqrt(Dxx.*Dxx+Dyy.*Dyy);
% px = Dxx./(mag+1e-10); py = Dyy./(mag+1e-10); 
% figure,quiver(py,px);

%%
%%%显示一阶导数
% Dx = imfilter(I,DGaussx,'conv');
% Dy = imfilter(I,DGaussy,'conv');
% mag = sqrt(Dx.*Dx+Dy.*Dy);
% px2 = Dx./(mag+1e-10); py2 = Dy./(mag+1e-10); 
% figure,quiver(py2,px2);
%%

%% using GVF Frangi Filter
% % [X,Y] = ndgrid(-round(3*Sigma):round(3*Sigma));
% % DGaussxx = 1/(2*pi*Sigma^4) * (X.^2/Sigma^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
% % DGaussxy = 1/(2*pi*Sigma^6) * (X .* Y)           .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
% % DGaussyy = DGaussxx';
% % Dxx = imfilter(I,DGaussxx,'conv');
% % Dxy = imfilter(I,DGaussxy,'conv');
% % Dyy = imfilter(I,DGaussyy,'conv');
% I = mapminmax(I, 0, 1);
% [Ex,Ey] = gradient(I);
% f = sqrt(Ex.*Ex+Ey.*Ey);
% [u,v] = GVF(f, 0.02, 100);
% % 
% % [line, row] = size(u);
% Y = [1    -1;
%      1    -1;
%      1    -1;
%        ];
% X = [-1 -1 -1;
%       1  1  1;
%            ];
% Dxx = imfilter(u,X,'conv');
% Dxy = imfilter(v,Y,'conv');
% Dyy = imfilter(u,Y,'conv');
% % for i = 2:line
% %     for j = 2:row
% %         Dxx(i,j) = u(i,j) - u(i-1,j);
% %         Dyy(i,j) = v(i,j) - v(i,j-1); 
% %         Dxy(i,j) = u(i,j) - u(i,j-1);
% %     end
% % end
% mag = sqrt(u.*u+v.*v);
% px2 = u./(mag+1e-10); py2 = v./(mag+1e-10); 
% figure,quiver(py2,px2);


function D = gradient(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
%
% Example:
%
% Fx = gradient3(F,'x');

[k,l] = size(F);
D  = zeros(size(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:) = (F(2,:) - F(1,:));
    D(k,:) = (F(k,:) - F(k-1,:));
    % Take centered differences on interior points
    D(2:k-1,:) = (F(3:k,:)-F(1:k-2,:))/2;
case 'y'
    D(:,1) = (F(:,2) - F(:,1));
    D(:,l) = (F(:,l) - F(:,l-1));
    D(:,2:l-1) = (F(:,3:l)-F(:,1:l-2))/2;
otherwise
    disp('Unknown option')
end

