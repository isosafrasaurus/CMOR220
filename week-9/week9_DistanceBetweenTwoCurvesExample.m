function ShapeAnalysisThurs2_Sol
load("DataClassification.mat")

for n = 1:100
    testdata(:, :, n) = preprocessData(testdata(:, :, n));
end
for n = 1:300
    trainingdata(:, :, n) = preprocessData(trainingdata(:, :, n));
end

Cnt=0;
for n=1:100
    TDist=zeros(300,100);
    TrClst=ceil(n/5);
    for m=1:300
        for k=1:100
            TstCv=[testdata(:,k:100,n),testdata(:,1:(k-1),n)];
            TDist(m,k)=DistanceCalc(TstCv,trainingdata(:,:,m));
        end
    end
    Val=min(min(TDist));
    [Idx,~]=find(TDist==Val);
    ObsClst=ceil(Idx/15);
    if TrClst==ObsClst
        Cnt=Cnt+1;
    end
end
CrctPcnt=(Cnt/100)*100;
display(CrctPcnt)
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

function preprocessedCurve = preprocessData(curve)
% Input: curve, a 2xNx1 matrix representing a curve
% Output: preprocessedCurve, a 2xNx1 matrix representing the preprocessed curve
mean_value = mean(curve, 2);  % Calculate the mean value of the curve
centeredCurve = curve - mean_value;  % Center the curve by subtracting the mean
preprocessedCurve = centeredCurve / norm(centeredCurve, 'fro');  % Normalize the centered curve
end