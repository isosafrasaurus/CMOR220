% Driver
function shapeAnalysis
% PART 1

% TRAINING SECTION
x = 0:0.02:pi;
y1 = [x;cos(x)];
y2 = [x;sin(x)];
y3 = [x;x.^2];
y4 = [x;x];

TrnData = zeros(2, length(x), 40);
for i=1:40
    if i<=10
        TrnData(:,:,i)=y1;
    elseif i>10 && i<=20
        TrnData(:,:,i)=y2;
    elseif i>20 && i<=30
        TrnData(:,:,i)=y3;
    else
        TrnData(:,:,i)=y4;
    end
end

for i=1:40
    TrnData(:,:,i)=GenNoise(TrnData(:,:,i));
    TrnData(:,:,i)=Preprocess(TrnData(:,:,i));
end

PWCM=zeros(40,40);

for n1=1:40
    for n2=1:40
        PWCM(n1,n2)=DistanceCalc(TrnData(:,:,n1),TrnData(:,:,n2));
    end
end

imagesc(PWCM); hold on;
colorbar;
title("pairwise distance matrix of the 40 open curves");

% TESTING SECTION

% create TstData, 25 per curve, 100 curves total
TstData = zeros(2, length(x), 100);
for i=1:100
    if i<=25
        TstData(:,:,i)=y1;
    elseif i>25 && i<=50
        TstData(:,:,i)=y2;
    elseif i>50 && i<=75
        TstData(:,:,i)=y3;
    else
        TstData(:,:,i)=y4;
    end
end

for i=1:100
    TstData(:,:,i)=GenNoise(TstData(:,:,i));
    TstData(:,:,i)=Preprocess(TstData(:,:,i));
end

%Classify Open Curves
Cnt=0;
for n1=1:100 % n1 refers to test data index
    TDist=zeros(1,40);
    TrClst=ceil(n1/25); % index of test data since they were made seqtly
    for n2=1:40 % use train data to classify
        TDist(n2)=DistanceCalc(TstData(:,:,n1),TrnData(:,:,n2));
    end
    [~,Idx]=min(TDist);
    ObsClst=ceil(Idx/10);
    if ObsClst==TrClst
        Cnt=Cnt+1;
    end
end
PercentCorrect=(Cnt/100)*100

clear;

% PART TWO

load("DataClassification.mat")
figure();
for n=1:20 %for each cluster
    subplot(4,5,n);
    for p=1:100
        plot(trainingdata(:,:,15*(n-1)+1));
    end
end

Cnt=0;
for n=1:5 % over test data
    TDist=zeros(40,200);
    TrClst=ceil(n/5);
    for m=1:15 % over train data
        for k=1:200
            TstCv=[TestData(:,k:200,n),TestData(:,1:(k-1),n)];
            TDist(m,k)=DistanceCalc(TstCv,TrainingData(:,:,m));
        end
    end
    Val=min(min(TDist));
    [Idx,~]=find(TDist==Val);
    ObsClst=ceil(Idx/10);
    if TrClst==ObsClst
        Cnt=Cnt+1;
    end
end
CrctPcnt=(Cnt/20)*100;
display(CrctPcnt)

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

function [D]=DistanceCalc(Cv1,Cv2)
A=Cv1*Cv2';
[U,~,V]=svd(A);
if det(U*V')>0
    D=norm(Cv1-(U*V')*Cv2,'fro');
else
    D=norm(Cv1-(U*[1 0; 0 -1]*V')*Cv2,'fro');
end
end