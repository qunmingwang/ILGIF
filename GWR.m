%%%%%GWR; Y and X are dependnt and independnt data in image format
%%%m: number of neighbors
%%%n:  number of variables
%%%x:  all independnt variabels, size m*(n+1)
%%%Wei:all weights, size m*m
%%%y:  all dependnt variabels, size m*1
%%%C: coefficients, size (n+1)*1

function [a,b]=GWR(X,Y,w);
m=(2*w+1)^2;
[Row,Col,n]=size(X);

for i=1:2*w+1
    for j=1:2*w+1
        D(i,j)=norm([i,j]-[w+1,w+1]);
    end
end
D_col=reshape(D,m,1);
Wei=zeros(m,m);
H=sqrt(2)*w;
for i=1:m
    Wei(i,i)=(1-(D_col(i)/H)^2)^2;
end

X_extend=Extend_cube(X,w);
Y_extend=Extend_plane(Y,w);
for i=w+1:Row+w
    for j=w+1:Col+w
        X_extend_local=X_extend(i-w:i+w,j-w:j+w,:);
        X_extend_local_D2=D3_D2(X_extend_local);
        x=[X_extend_local_D2',ones(m,1)];
        Y_extend_local=Y_extend(i-w:i+w,j-w:j+w);
        y=reshape(Y_extend_local,m,1);
        C(i-w,j-w,:)=inv(x'*Wei*x)*x'*Wei*y;
    end
end
a=C(:,:,1:n);
b=C(:,:,n+1);