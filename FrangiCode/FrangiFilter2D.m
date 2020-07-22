function [outIm,Direction, whatScale] = FrangiFilter2D(I, sigmas, options)
defaultoptions = struct('FrangiScaleRange', [1 10], 'FrangiScaleRatio', 2, 'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 10, 'verbose',true,'BlackWhite',true);%0.5 15

if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('FrangiFilter2D:unknownoption','unknown options found');
    end
end
%sigmas = [3];
beta  = 2*options.FrangiBetaOne^2;
c     = 2*options.FrangiBetaTwo^2;
ALLfiltered=zeros([size(I) length(sigmas)]);
ALLangles=zeros([size(I) length(sigmas)]);
[line, row] = size(I);
for i = 1:length(sigmas),
    if(options.verbose)
        
    end
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));
    Dxx = (sigmas(i)^2)*Dxx;
    Dxy = (sigmas(i)^2)*Dxy;
    Dyy = (sigmas(i)^2)*Dyy;
    [Lambda2,Lambda1,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);
%     angles =atan2(Ix,Iy);
    angles = zeros(line,row);
    for ii = 1:line
        for jj = 1:row
            angles(ii, jj) = atan(-((Iy(ii,jj)/Ix(ii,jj))));
        end
    end
    Lambda1(Lambda1==0) = eps;
    Rb = (Lambda2./Lambda1).^2;%
    S2 =  Lambda1.^2+Lambda2.^2;% 

   Ifiltered =  exp(-Rb/beta).*(ones(size(I))-exp(-S2/c));
   
   if(options.BlackWhite)
        Ifiltered(Lambda1<0)=0;
    else
        Ifiltered(Lambda1>0)=0;
    end
    ALLfiltered(:,:,i) = Ifiltered;
    ALLangles(:,:,i) = angles;
end

if length(sigmas) > 1,
    [outIm,whatScale] = max(ALLfiltered,[],3);
    outIm = reshape(outIm,size(I));
    if(nargout>1)
        whatScale = reshape(whatScale,size(I));
    end
    %if(nargout>2)
        Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
    %end
else
    outIm = reshape(ALLfiltered,size(I));
    if(nargout>1)
            whatScale = ones(size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles,size(I));
    end
end
