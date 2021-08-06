
function rba2s2fModel(nsamples)
%This function find the best fitting 2 slope model for RBA data and the
%uncertainty on each slope and breakpoint.
%The two slopes of the model are free (instead of having one fixed slope)

%nsamples is the number of samples used for calculating the 
%uncertainty of the breakpoint and slope with a bootstrapping 
%method.
%% load data

%RBA data
newRBA=readtable('RBArobbins18Goossens20_rho2550_taper.csv');

%sort RBA by diameter 
[~,indxSort]=sort(newRBA.Diam);
newRBA=newRBA(indxSort,:);

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

% plot cumulative probability of breakpoint and slopes

figure(1)
subplot(3,1,1)
plot(xbreakValues,cumprobBreak)
xlabel('Breakpoint (km)')
ylabel('Cumulative probability')
hold on
line([min(allBreakPoints) max(allBreakPoints)],[0.025 0.025],'LineStyle','--',...
      'Color','k')
hold on
line([min(allBreakPoints) max(allBreakPoints)],[0.975 0.975],'LineStyle','--',...
      'Color','k')
xlim([min(allBreakPoints) max(allBreakPoints)])
title(['Breakpoint=' num2str(bestBreak) ' +'...
       num2str(upperBreakpoint-bestBreak) ' -' ...
       num2str(bestBreak-lowerBreakpoint)])
   
subplot(3,1,2)
plot(xS1Values,cumprobS1)
xlabel('Slope 1 (mGal/km)')
ylabel('Cumulative probability')
hold on
line([min(ibestCoeff1(:,1)) max(ibestCoeff1(:,1))],[0.025 0.025],'LineStyle','--',...
      'Color','k')
hold on
line([min(ibestCoeff1(:,1)) max(ibestCoeff1(:,1))],[0.975 0.975],'LineStyle','--',...
      'Color','k')
xlim([min(ibestCoeff1(:,1)) max(ibestCoeff1(:,1))])
title(['Slope 1=' num2str(bestCoeff1(1)) ' +'...
       num2str(upperS1-bestCoeff1(1)) ' -' ...
       num2str(bestCoeff1(1)-lowerS1)])

subplot(3,1,3)
plot(xS2Values,cumprobS2)
xlabel('Slope 2 (mGal/km)')
ylabel('Cumulative probability')
hold on
line([min(ibestCoeff2(:,1)) max(ibestCoeff2(:,1))],[0.025 0.025],'LineStyle','--',...
      'Color','k')
hold on
line([min(ibestCoeff2(:,1)) max(ibestCoeff2(:,1))],[0.975 0.975],'LineStyle','--',...
      'Color','k')
xlim([min(ibestCoeff2(:,1)) max(ibestCoeff2(:,1))])
title(['Slope 2=' num2str(bestCoeff2(1)) ' +'...
       num2str(upperS2-bestCoeff2(1)) ' -' ...
       num2str(bestCoeff2(1)-lowerS2)])

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


%% plot data   

%make coordinates for fill in area Slopes
yCoorLimitsF=[lowerF fliplr(upperF)];
xCoorLimitsF=[newRBA.Diam' fliplr(newRBA.Diam')];


figure(2)

subplot(2,1,1)

limRBA=100;
plot(newRBA.Diam,newRBA.RBA,'.')
ylim([-limRBA limRBA])
xlabel('Crater diameter (km)')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorLimitsF,yCoorLimitsF,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(newRBA.Diam,bestF,...
    'Color','k',...
    'LineWidth',3)
hold on
line([bestBreak bestBreak],...
     [-limRBA,limRBA],'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
hold on
line([lowerBreakpoint lowerBreakpoint],...
     [-limRBA,limRBA],...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
hold on
line([upperBreakpoint upperBreakpoint],...
     [-limRBA,limRBA],...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)


%zoom in-------------------
subplot(2,1,2)

plot(newRBA.Diam,newRBA.RBA,'.')
ylim([-10 10])
xlabel('Crater diameter (km)')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorLimitsF,yCoorLimitsF,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(newRBA.Diam,bestF,...
    'Color','k',...
    'LineWidth',3)
hold on
line([bestBreak bestBreak],...
     [-limRBA,limRBA],'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
hold on
line([lowerBreakpoint lowerBreakpoint],...
     [-limRBA,limRBA],...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
hold on
line([upperBreakpoint upperBreakpoint],...
     [-limRBA,limRBA],...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)



end
