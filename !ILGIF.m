%%%%This is the code for ILGIF produced by Dr Qunming Wang; Email: wqm11111@126.com
%%%%Copyright belong to Qunming Wang
%%%%When using the code, please cite the fowllowing paper
%%%%Q. Wang, W. Shi, P. M. Atkinson. Information loss-guided multi-resolution image fusion. IEEE Transactions on Geoscience and Remote Sensing, 2020, 58(1): 45¨C57.

clear all;
load S2_20m;%%%20m Sentinel-2 bands in a image cube (6 bands of reflectance)
load S2_10m;%%%10m Sentinel-2 bands in a image cube (4 bands of reflectance)
s=2;
I_MS=S2_20m;
I_PAN=S2_10m;

w=1;
sigma=s/2;
PSFh=PSF_template(s,w,sigma);%%%Gaussian PSF; sigma was set to 0.5*s here, but needs to be specified in practice

Sill_min=1;
Range_min=0.5;
L_sill=20;
L_range=20;
rate=0.1;
H=20;

I_PAN_L=dowmsample_cube(I_PAN,s,w,PSFh);%%%upscaling the fine (e.g., 10 m) bands to coarse ones (e.g., 20 m)

for i=1:size(I_PAN,3)
    PAN_DS(:,:,i)=ATPK_DS(I_PAN_L(:,:,i),s,Sill_min,Range_min,L_sill,L_range,rate,H,w,PSFh);%%%downscaling the simulated coarse bands
end

ILoriginal_cube=I_PAN-PAN_DS;%%%Information loss of the fine bands

w2=3;

tic
for i=1:size(I_MS,3)
    MS_DS(:,:,i)=ATPK_DS(I_MS(:,:,i),s,Sill_min,Range_min,L_sill,L_range,rate,H,w,PSFh);%%%downscaling the observed coarse bands
    %%%GWR
    [a,b]=GWR(I_PAN_L,I_MS(:,:,i),w2);%%%coefficients estimated by GWR
    IL=sum(imresize(a,s,'nearest').*ILoriginal_cube,3);%%%Information loss estimated for the observed coarse bands
    Z(:,:,i)=MS_DS(:,:,i)+IL;
end
alltime=toc

FalseColorf=Z(:,:,[3,2,1]);xf=imadjust(FalseColorf,stretchlim(FalseColorf),[]);figure,imshow(xf);
