function I=imgaussian(I,sigma,siz)
if(~exist('siz','var')), siz=sigma*6; end
if(sigma>0)
    x=-ceil(siz/2):ceil(siz/2);
    H = exp(-(x.^2/(2*sigma^2)));
    H = H/sum(H(:));

    if(ndims(I)==1)
        I=imfilter(I,H, 'same' ,'replicate');
    elseif(ndims(I)==2)
        Hx=reshape(H,[length(H) 1]);
        Hy=reshape(H,[1 length(H)]);
        I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
    elseif(ndims(I)==3)
        if(size(I,3)<4) % Detect if 3D or color image
            Hx=reshape(H,[length(H) 1]);
            Hy=reshape(H,[1 length(H)]);
            for k=1:size(I,3)
                I(:,:,k)=imfilter(imfilter(I(:,:,k),Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
            end
        else
            Hx=reshape(H,[length(H) 1 1]);
            Hy=reshape(H,[1 length(H) 1]);
            Hz=reshape(H,[1 1 length(H)]);
            I=imfilter(imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
        end
    else
        error('imgaussian:input','unsupported input dimension');
    end
end