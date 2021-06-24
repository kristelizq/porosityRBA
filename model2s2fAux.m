
function [bestCoeff1,bestCoeff2,bestBreak,bestF,...
          upperBreakpoint,lowerBreakpoint,...
          upperS1,lowerS1,upperS2,lowerS2,...
          xCoorLimitsF,yCoorLimitsF]=model2s2fAux(rbaData,nsamples)
%This function find the best fitting 2 slope model for RBA data and the
%uncertainty on each slope and breakpoint.
%The two slopes of the model are free (instead of having one fixed slope)

%% load data

%RBA data
newRBA=rbaData;

%% Find best fit model

allBreakPoints=round(min(newRBA.Diam)):0.2:round(max(newRBA.Diam));
nbreaks=numel(allBreakPoints);

n=numel(newRBA.Diam);
p=5;
BIC=NaN(nbreaks,1);
coef_all_1=zeros(nbreaks,2);
coef_all_2=zeros(nbreaks,2);
f_all=zeros(nbreaks,n);

for i=3:nbreaks-3
    
    %find data points before break
    indxDiamLessBreak=find(newRBA.Diam<=allBreakPoints(i));
    coef_all_1(i,:) = polyfit(newRBA.Diam(indxDiamLessBreak),...
                      newRBA.RBA(indxDiamLessBreak),1);
    f_all(i,indxDiamLessBreak)=coef_all_1(i,1)*...
                            newRBA.Diam(indxDiamLessBreak)+...
                            coef_all_1(i,2);
                        
    %get slope of second line constrained to matching one data point
    % in f_2s2f

    V=[newRBA.Diam(indxDiamLessBreak(end):n),...
       ones(numel(newRBA.Diam(indxDiamLessBreak(end):n)),1)];
    inx=[];
    iny=[];
    ix=newRBA.Diam(indxDiamLessBreak(end)).^(1:-1:0);
    iy=f_all(i,indxDiamLessBreak(end));
    
    coef_all_2(i,:)=lsqlin(V,f_all(i,indxDiamLessBreak(end):end),...
                              inx,iny,ix,iy);

    %evaluate coefficients obtained
    f_all(i,indxDiamLessBreak(end):n)=coef_all_2(i,1)*...
                            newRBA.Diam(indxDiamLessBreak(end):n)+...
                            coef_all_2(i,2);
    res=(f_all(i,:)-newRBA.RBA')*(f_all(i,:)-newRBA.RBA')';
    L=-n/2*log(res);
    BIC(i)=L-1/2*p*log(n);
    
end

[maxBIC,indxMaxBIC]=max(BIC);
 
bestBreak=allBreakPoints(indxMaxBIC);
bestCoeff1=coef_all_1(indxMaxBIC,:);
bestCoeff2=coef_all_2(indxMaxBIC,:);
bestF=f_all(indxMaxBIC,:);

%% get uncertainty in breakpoint and slopes


ibestBreak=zeros(nsamples,1);
ibestCoeff1=zeros(nsamples,2);
ibestCoeff2=zeros(nsamples,2);
ibestFsub=zeros(nsamples,n-1);

for isample=1:nsamples
    
    RBAsample=datasample(newRBA,n-1,1);
    [~,I]=sort(RBAsample.Diam);
    RBAsample=RBAsample(I,:);
    
    BIC=NaN(nbreaks,1);
    coef_sub_1=zeros(nbreaks,2);
    coef_sub_2=zeros(nbreaks,2);
    f_sub=zeros(nbreaks,n-1);

    for i=3:nbreaks-3

        %find data points before break
        indxDiamLessBreak=find(RBAsample.Diam<=allBreakPoints(i));
        coef_sub_1(i,:) = polyfit(RBAsample.Diam(indxDiamLessBreak),...
                          RBAsample.RBA(indxDiamLessBreak),1);
        f_sub(i,indxDiamLessBreak)=coef_sub_1(i,1)*...
                                RBAsample.Diam(indxDiamLessBreak)+...
                                coef_sub_1(i,2);

        %get slope of second line constrained to matching one data point
        % in f_2s2f

        V=[RBAsample.Diam(indxDiamLessBreak(end):n-1),...
           ones(numel(RBAsample.Diam(indxDiamLessBreak(end):n-1)),1)];
        inx=[];
        iny=[];
        ix=RBAsample.Diam(indxDiamLessBreak(end)).^(1:-1:0);
        iy=f_sub(i,indxDiamLessBreak(end));

        coef_sub_2(i,:)=lsqlin(V,f_sub(i,indxDiamLessBreak(end):end),...
                                  inx,iny,ix,iy);

        %evaluate coefficients obtained
        f_sub(i,indxDiamLessBreak(end):n-1)=coef_sub_2(i,1)*...
                                RBAsample.Diam(indxDiamLessBreak(end):n-1)+...
                                coef_sub_2(i,2);
        res=(f_sub(i,:)-RBAsample.RBA')*(f_sub(i,:)-RBAsample.RBA')';
        L=-(n-1)/2*log(res);
        BIC(i)=L-1/2*p*log(n-1);

    end

    [maxBIC,indxMaxBIC]=max(BIC);
    ibestBreak(isample)=allBreakPoints(indxMaxBIC);
    ibestCoeff1(isample,:)=coef_sub_1(indxMaxBIC,:);
    ibestCoeff2(isample,:)=coef_sub_2(indxMaxBIC,:);
    ibestFsub(isample,:)=f_sub(indxMaxBIC,:);
    isample
end

%get limits breakpoint and slopes----------------------------
%breakpoint
[cumprobBreak,xbreakValues] = ecdf(ibestBreak);

%find lower limit breakpoint
indxlowerBreak=max(find(cumprobBreak<=0.025));
lowerBreakpoint=xbreakValues(indxlowerBreak);
%find upper limit breakpoint
indxupperBreak=min(find(cumprobBreak>=0.975));
upperBreakpoint=xbreakValues(indxupperBreak);

%slope 1
[cumprobS1,xS1Values] = ecdf(ibestCoeff1(:,1));

%find lower limit slope 1
indxlowerS1=max(find(cumprobS1<=0.025));
lowerS1=xS1Values(indxlowerS1);
%find upper limit slope 1
indxupperS1=min(find(cumprobS1>=0.975));
upperS1=xS1Values(indxupperS1);

% slope 2
[cumprobS2,xS2Values] = ecdf(ibestCoeff2(:,1));

%find lower limit slope 2
indxlowerS2=max(find(cumprobS2<=0.025));
lowerS2=xS2Values(indxlowerS2);
%find upper limit slope 2
indxupperS2=min(find(cumprobS2>=0.975));
upperS2=xS2Values(indxupperS2);

%% get limits of F function (slopes and interecepts)

lowerF=zeros(1,n);
upperF=zeros(1,n);

for i=1:n
    
    [cumprobF,xFvalues] = ecdf(f_all(:,i));

    %find lower limit slope 2
    indxlowerF=max(find(cumprobF<=0.025));
    lowerF(i)=xFvalues(indxlowerF);
    %find upper limit slope 2
    indxupperF=min(find(cumprobF>=0.975));
    upperF(i)=xFvalues(indxupperF);
    
end
  

%make coordinates for fill in area Slopes
yCoorLimitsF=[lowerF fliplr(upperF)];
xCoorLimitsF=[newRBA.Diam' fliplr(newRBA.Diam')];


end