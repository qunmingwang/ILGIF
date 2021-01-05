function [RMSE,CC,ERGAS,UIQI,SAM]=evaluate_relation_spectral(realdata,predictdata,s);
[a,b,c]=size(realdata);
for i=1:c
    P=predictdata(:,:,i);R=realdata(:,:,i);
    
    RMSE1=sum(sum((R-P).^2));RMSE0(i)=sqrt(RMSE1/(a*b));

    C_1=sum(sum(P.*R))-a*b*mean(mean(P))*mean(mean(R));C_2=sum(sum(P.^2))-a*b*mean(mean(P))^2;C_3=sum(sum(R.^2))-a*b*mean(mean(R))^2;
    CC0(i)=C_1/sqrt(C_2*C_3);
    
    ERGAS0(i)=RMSE0(i)/mean(mean(R));
    
    UIQI_1=4*mean(mean(P))*mean(mean(R))*C_1;%%%/(a*b)
    UIQI_2=(mean(mean(P))^2+mean(mean(R))^2)*(C_2+C_3);%%%/(a*b)
    UIQI0(i)=UIQI_1/UIQI_2;
end
RMSE=[RMSE0,mean(RMSE0)];
CC=[CC0,mean(CC0)];
ERGAS=100*norm(ERGAS0)/(s*sqrt(c));
UIQI=[UIQI0,mean(UIQI0)];

for i=1:a
    for j=1:b
        VP=D3_D2(predictdata(i,j,:));VR=D3_D2(realdata(i,j,:));
        
        SAM0(i,j)=VP'*VR/((norm(VP)+0.0000001)*(norm(VR)+0.0000001));
        
    end
end
SAM=mean(mean(acos(SAM0)));