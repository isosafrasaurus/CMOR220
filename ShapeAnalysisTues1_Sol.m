function ShapeAnalysisTues1_Sol
x=0:0.01:4;
y1=[x;7*sin(x)];
y2=[x;0.5*exp(x-2)]
y3=[x;4*sqrt(x+5)]
y4=[x;0.5*log(x+6)];
TrnData=zeros(2,length(x),100);
for n=1:20
    TrnData(:,:,n)=GenNoise(y1);
    TrnData(:,:,n+20)=GenNoise(y2);
    TrnData(:,:,n+40)=GenNoise(y3);
    TrnData(:,:,n+60)=GenNoise(y4);
end
for n=1:100
    TrnData(:,:,n)=Preprocess(TrnData(:,:,n));
end
PWCM=zeros(100,100); %Preallocate PWCM
for n=1:100
    for m=1:100
        PWCM(n,m)=DistanceCalculator(TrnData(:,:,n),TrnData(:,:,m));
    end
end
%Need to image the matrix here
figure
imagesc(PWCM) ; colorbar ; hold on ; grid on
title('\fontname{TW Cen MT} PWCM Matrix of the Open Curves')
end

function [Distance]=DistanceCalculator(Cv1,Cv2)
A=Cv1*Cv2';
[U,~,V]=svd(A);
if det(U*V')>0
    Distance=norm(Cv1-(U*V')*Cv2,'fro');
else
    Distance=norm(Cv1-(U*[1 0; 0 -1]*V')*Cv2,'fro');
end
end

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