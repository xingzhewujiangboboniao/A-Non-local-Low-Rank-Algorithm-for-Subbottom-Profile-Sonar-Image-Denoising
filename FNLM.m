function  [Output]=FNLM(x,f,t,h, Ivessel, l, s)
warning off;

oldx = x;
x = Ivessel;

[m n]=size(x);
% Memory for the output
Output=[];

% Replicate the boundaries of the x image

D_ori = padarray(x,[f+t f+t],'symmetric');
oldx = padarray(oldx,[f+t f+t],'symmetric');
kernel=make_kernel(f);
gsqu_sum2 = conv2(oldx.*oldx,kernel,'same');
h=h^2;

parfor i=1:m
    i
    for j=1:n
        i1 = i+ t+f;
        j1 = j+ t+f;
        
        rw2= oldx(i1-f:i1+f , j1-f:j1+f); % reference window
        rw2(:,1:end)=rw2(:,2*f+1:-1:1);
        rw2(1:end,:)=rw2(2*f+1:-1:1,:);
        rw2=rw2.*kernel;
        
        bw2 = oldx(i1-t-f:i1+t+f,j1-t-f:j1+t+f); 
        cv_bw2 = conv2(bw2,rw2,'valid');
        %----------------------------------------------------------
        
        gsq_dis2=gsqu_sum2(i1,j1)+gsqu_sum2(i1-t:i1+t,j1-t:j1+t)-2*cv_bw2;  % Eqn. (5) in the paper;
        weight2=exp(-gsq_dis2/h);
        weight2(t+1,t+1)=0;
        weight2(t+1,t+1)=max(weight2(:));
        weight2 = weight2/max(weight2(:));
        %%
        weight2(:,t+1) = 0;

        wei = zeros(size(weight2));
        Dir = D_ori(i1-t:i1+t,j1-t:j1+t);
        for ii = 1:size(weight2,1)
            for jj = 1:size(weight2,2)
                a = ii - (t+1);
                b = jj - (t+1);
                if (ii~=t+1 || jj~=t+1)
                    Pont_angle = atan(a/b);
                    E = abs(Dir(t+1,t+1) - Pont_angle);%when the centre direction is wrong, the result is bad.
                    lamada = 0.3;
                    if (E> lamada)
                        wei(ii,jj ) = 0;%0.0001;
                        continue;
                    end
                    wei(ii,jj) = (1-(E/lamada).^2).^2 ;%exp(-abs(wei(ii,jj) - Dir(ii,jj)));
                end
            end
        end
        wei(t+1,t+1) = 1;
        wei2 = zeros(size(weight2));
        for ii = 1:size(weight2,1)
            for jj = 1:size(weight2,2)
                a = (t+1) - ii;
                b = (t+1) - jj;
                if (ii~=t+1 || jj~=t+1)
                    Pont_angle2 = atan(a/b);
                    E2 = abs(Dir(ii,jj) - Pont_angle2);
                    lamada = 0.3;
                    if (E2> lamada)
                        wei2(ii,jj ) = 0;%0.0001;
                        continue;
                    end
                        wei2(ii,jj) = (1-(E2/lamada).^2).^2 ;%exp(-abs(wei(ii,jj) - Dir(ii,jj)));
                end
            end
        end

        W =  wei.*wei2.*weight2;%   wei.*  abs(wei-ones(size(wei))).*   weight.*D_2.*wei.*.*weight2
        if (sum(sum(wei2.*weight2)) < 2*sum(sum(wei.*wei2.*weight2)))
            W =  wei.*wei2.*weight2;
        else
            W =  wei2.*weight2;
        end
        [ L,~ ] =lowrank(W,oldx, i1, j1, t, f, l, s); 
        Output(i,j) = mean(L(floor(size(L,1)/2)+1,:));
    end
end

function [kernel] = make_kernel(f)

kernel=zeros(2*f+1,2*f+1);
for d=1:f
    value= 1 / (2*d+1)^2 ;
    for i=-d:d
        for j=-d:d
            kernel(f+1-i,f+1-j)= kernel(f+1-i,f+1-j) + value ;
        end
    end
end
kernel = kernel ./ f;



