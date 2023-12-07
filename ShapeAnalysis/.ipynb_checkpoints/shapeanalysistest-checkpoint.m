function shapeanalysistest
y1=[x;7*sin(x)];
y2=[x;0.5*exp(x-2)]
y3=[x;4*sqrt(x+5)]
y4=[x;0.5*log(x+6)];

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