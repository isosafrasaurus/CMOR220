function ShapeAnalysisThurs1_Sol
t=0:0.01:5;
CubR=[t; t.^(1/3)];
E1=[t; exp(t)];
E2=[t; exp(t.^2)];
Sine2=[t; sin(2*t)];
X=[t;t];
TrnData=zeros(2,length(t),50);
for i=1:50
if i<=10
TrnData(:,:,i)=CubR;
elseif i>10 && i<=20
TrnData(:,:,i)=E1;
elseif i>20 && i<=30
TrnData(:,:,i)=E2;
elseif i>30 && i<=40
TrnData(:,:,i)=Sine2;
else
TrnData(:,:,i)=X;
end
end
%--------------------------------------------------------------------------
%DATA SETUP AND GENERATION. PRACTICE PROBLEM BEGINS BELOW
%DRIVER
for i=1:50
TrnData(:,:,i)=GenNoise(TrnData(:,:,i));
TrnData(:,:,i)=Preprocess(TrnData(:,:,i));
end
PWCM=zeros(50,50);
for n1=1:50
for n2=1:50
PWCM(n1,n2)=DistanceCalc(TrnData(:,:,n1),TrnData(:,:,n2));
end
end
figure
imagesc(PWCM) ; hold on ;
colorbar;
title("\fontname{TW Cen MT} ANOTHER PWCM: MORE CURVES, MORE VALUE")
subtitle("\fontname{TW Cen MT} PRICE: $2999.99",'FontWeight','bold')
end
%FUNCTIONS
function [nC]=GenNoise(C)
for n=1:length(C(1,:,:))
nC(1,n,:)=C(1,n,:)+randn()*0.05;
nC(2,n,:)=C(2,n,:)+randn()*0.05;
end
end

function [pC]=Preprocess(C)
[M]=mean(C,2);
centerC=C-M;
pC=centerC/norm(centerC,'fro');
end

function [D]=DistanceCalc(Cv1,Cv2)
A=Cv1*Cv2';
[U,~,V]=svd(A);
if det(U*V')>0
D=norm(Cv1-(U*V')*Cv2,'fro');
else
D=norm(Cv1-(U*[1 0; 0 -1]*V')*Cv2,'fro');
end
end