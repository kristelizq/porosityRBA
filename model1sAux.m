
function [bestCoeff1,bestF,...
          upperS1,lowerS1,...
          xCoorLimitsF,yCoorLimitsF]=model1sAux(rbaData,nsamples)
%This function find the best fitting 1 slope model for RBA data and the
%uncertainty on the slope.

%% load data

%RBA data
newRBA=rbaData;
n=numel(newRBA.Diam);

%% Find best fit model

%get coefficients
bestCoeff1 = polyfit(newRBA.Diam,newRBA.RBA,1);
bestF=bestCoeff1(1)*newRBA.Diam+bestCoeff1(2);

%% get uncertainty in slope

ibestCoeff1=zeros(nsamples,2);
ibestF=zeros(nsamples,n-1);

for isample=1:nsamples
    
    RBAsample=datasample(newRBA,n-1,1);
    [~,I]=sort(RBAsample.Diam);
    RBAsample=RBAsample(I,:);
    
    ibestCoeff1(isample,:) = polyfit(RBAsample.Diam,RBAsample.RBA,1);
    ibestF(isample,:)=ibestCoeff1(isample,1)*...
                      RBAsample.Diam+ibestCoeff1(isample,2);
    
    isample
end

%get limits slope----------------------------

[cumprobS1,xS1Values] = ecdf(ibestCoeff1(:,1));

%find lower limit slope 1
indxlowerS1=max(find(cumprobS1<=0.025));
lowerS1=xS1Values(indxlowerS1);
%find upper limit slope 1
indxupperS1=min(find(cumprobS1>=0.975));
upperS1=xS1Values(indxupperS1);


%% get limits of F function (slopes and interecepts)

lowerF=zeros(1,n-1);
upperF=zeros(1,n-1);

for i=1:n-1
    
    [cumprobF,xFvalues] = ecdf(ibestF(:,i));

    %find lower limit slope 2
    indxlowerF=max(find(cumprobF<=0.025));
    lowerF(i)=xFvalues(indxlowerF);
    %find upper limit slope 2
    indxupperF=min(find(cumprobF>=0.975));
    upperF(i)=xFvalues(indxupperF);
    
end
  

%make coordinates for fill in area Slopes
yCoorLimitsF=[lowerF fliplr(upperF)];
xCoorLimitsF=[newRBA.Diam(1:end-1)' fliplr(newRBA.Diam(1:end-1)')];


end