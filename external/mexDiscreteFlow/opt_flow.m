function [vx,vy,energylist] = opt_flow(im1,im2,params,isdisplay,Segmentation)

if nargin < 3
    alpha=0.01; % smooth coeff
    d=alpha*20; % smooth thres
    eta=0.001; % small displacement coeff
    nlevels=4; % n levels of pyramid
    wsize=3; % search window size below top
    topwsize=10; % search window size of top
    nIterations=40; % n iters below top
    nTopIterations=100; % n iters of top
else
    alpha=params.alpha;
    d=params.d;
    eta=params.eta;
    nlevels=params.nlevels;
    wsize=params.wsize;
    topwsize=params.topwsize;
    nIterations=params.nIterations;
    nTopIterations=params.nTopIterations;
end
    
if exist('isdisplay','var')~=1
    isdisplay=false;
end

if exist('Segmentation','var')==1
    IsSegmentation=true;
else
    IsSegmentation=false;
end

% build the pyramid
pyrd(1).im1=im1;
pyrd(1).im2=im2;
if IsSegmentation
    pyrd(1).seg=Segmentation;
end

for i=2:nlevels
    pyrd(i).im1=imresize(imfilter(pyrd(i-1).im1,fspecial('gaussian',5,0.67),'same','replicate'),0.5,'bicubic');
    pyrd(i).im2=imresize(imfilter(pyrd(i-1).im2,fspecial('gaussian',5,0.67),'same','replicate'),0.5,'bicubic');
    if IsSegmentation
        pyrd(i).seg=imresize(pyrd(i-1).seg,0.5,'nearest');
    end
end

for i=1:nlevels
    [height,width,nchannels]=size(pyrd(i).im1);
    [height2,width2,nchannels]=size(pyrd(i).im2);
    [xx,yy]=meshgrid(1:width,1:height);    
    pyrd(i).xx=round((xx-1)*(width2-1)/(width-1)+1-xx);
    pyrd(i).yy=round((yy-1)*(height2-1)/(height-1)+1-yy);
end

nIterationArray=round(linspace(nIterations,nIterations,nlevels));

% main loop
for i=nlevels:-1:1
    if isdisplay
        fprintf('Level: %d...',i);
    end
    [height,width,nchannels]=size(pyrd(i).im1);
    [height2,width2,nchannels]=size(pyrd(i).im2);
    [xx,yy]=meshgrid(1:width,1:height);
    
    if i==nlevels % top
        vx=pyrd(i).xx;
        vy=pyrd(i).yy;
        
        winSizeX=ones(height,width)*topwsize;
        winSizeY=ones(height,width)*topwsize;
    else
        vx=round(pyrd(i).xx+imresize(vx-pyrd(i+1).xx,[height,width],'bicubic')*2);
        vy=round(pyrd(i).yy+imresize(vy-pyrd(i+1).yy,[height,width],'bicubic')*2);
        
        winSizeX=ones(height,width)*(wsize+i-1);
        winSizeY=ones(height,width)*(wsize+i-1);
    end

    Im1=pyrd(i).im1;
    Im2=pyrd(i).im2;

    % compute the image-based coefficient
    if IsSegmentation
        imdiff=zeros(height,width,2);
        imdiff(:,1:end-1,1)=double(pyrd(i).seg(:,1:end-1)==pyrd(i).seg(:,2:end));
        imdiff(1:end-1,:,2)=double(pyrd(i).seg(1:end-1,:)==pyrd(i).seg(2:end,:));
        Im_s=imdiff*alpha+(1-imdiff)*alpha*0.01;
        Im_d=imdiff*alpha*100+(1-imdiff)*alpha*0.01*20;
    end

    if i==nlevels
        if IsSegmentation
            [flow,foo]=mexDiscreteFlow(Im1,Im2,[alpha,d,eta*2^(i-1),nTopIterations,2,topwsize],vx,vy,winSizeX,winSizeY,Im_s,Im_d);
        else
            [flow,foo]=mexDiscreteFlow(Im1,Im2,[alpha,d,eta*2^(i-1),nTopIterations,2,topwsize],vx,vy,winSizeX,winSizeY);
        end
    else
        if IsSegmentation
            [flow,foo]=mexDiscreteFlow(Im1,Im2,[alpha,d,eta*2^(i-1),nIterationArray(i),nlevels-i,wsize],vx,vy,winSizeX,winSizeY,Im_s,Im_d);
        else
            [flow,foo]=mexDiscreteFlow(Im1,Im2,[alpha,d,eta*2^(i-1),nIterationArray(i),nlevels-i,wsize],vx,vy,winSizeX,winSizeY);
        end
    end
    energylist(i).data=foo;
    vx=flow(:,:,1);
    vy=flow(:,:,2);
    if isdisplay
        fprintf('done!\n');
    end
end


function winSizeIm=decideWinSize(offset,wsize)

% design the DOG filter
f1=fspecial('gaussian',9,1);
f2=fspecial('gaussian',9,.5);
f=f2-f1;

foo=imfilter(abs(imfilter(offset,f,'same','replicate')),fspecial('gaussian',9,1.5),'same','replicate');

Min=min(foo(:));
Max=max(foo(:));
winSizeIm=wsize*(foo-Min)/(Max-Min)+wsize;


