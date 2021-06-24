function BICmodels()
%This function calculates the Bayesian Information Criterion (BIC)
%of models with 1 slope, 2 free slopes and 2 slopes with 1 fixed slope.

newRBA=readtable('RBArobbins18Goossens20_taper.csv');
%sort RBA by diameter 
[~,indxSort]=sort(newRBA.Diam);
newRBA=newRBA(indxSort,:);

n=numel(newRBA.Diam);
%% one slope model

%get coefficients
coefFit_1s = polyfit(newRBA.Diam,newRBA.RBA,1);
f_1s=coefFit_1s(1)*newRBA.Diam+coefFit_1s(2);

%calculate maximized likelihood L
p_1s=3; %number of degrees of freedom
res_1s=(f_1s-newRBA.RBA)'*(f_1s-newRBA.RBA); %residual
L_1s=-n/2*log(res_1s); %likelihood
BIC_1s=L_1s-1/2*p_1s*log(n); % Bayesian Information Criterion



%% two slope model (1 fixed)
%do this for each breakpoint value
breakp=round(min(newRBA.Diam)):1:round(max(newRBA.Diam));
nbreaks=numel(breakp);

p_2s=4;
BIC_2s=NaN(nbreaks,1);
coefFit_2sAll=zeros(nbreaks,2);
f_2sAll=zeros(nbreaks,n);

for i=3:nbreaks-3
    
    %find data points before break
    indxDiamLessBreak=find(newRBA.Diam<=breakp(i));
    coefFit_2sAll(i,:) = polyfit(newRBA.Diam(indxDiamLessBreak),...
                      newRBA.RBA(indxDiamLessBreak),1);
    f_2sAll(i,indxDiamLessBreak)=coefFit_2sAll(i,1)*...
                            newRBA.Diam(indxDiamLessBreak)+...
                            coefFit_2sAll(i,2);
    f_2sAll(i,indxDiamLessBreak(end)+1:n)=polyval(coefFit_2sAll(i,:),...
                                  breakp(i));
    res_2s=(f_2sAll(i,:)-newRBA.RBA')*(f_2sAll(i,:)-newRBA.RBA')';
    L_2s=-n/2*log(res_2s);
    BIC_2s(i)=L_2s-1/2*p_2s*log(n);
    
end

[~,indxmaxBIC2]=max(BIC_2s);

breakpoint_2s=breakp(indxmaxBIC2);


%% two slope model 2 free

p_2s2f=5;
BIC_2s2f=NaN(nbreaks,1);
coefFit_2s2f_1All=zeros(nbreaks,2);
coefFit_2s2f_2All=zeros(nbreaks,2);
f_2s2fAll=zeros(nbreaks,n);

for i=3:nbreaks-3
    
    %find data points before break
    indxDiamLessBreak=find(newRBA.Diam<=breakp(i));
    coefFit_2s2f_1All(i,:) = polyfit(newRBA.Diam(indxDiamLessBreak),...
                      newRBA.RBA(indxDiamLessBreak),1);
    f_2s2fAll(i,indxDiamLessBreak)=coefFit_2s2f_1All(i,1)*...
                            newRBA.Diam(indxDiamLessBreak)+...
                            coefFit_2s2f_1All(i,2);
                        
    %get slope of second line constrained to matching one data point
    % in f_2s2f

    V=[newRBA.Diam(indxDiamLessBreak(end):n),...
       ones(numel(newRBA.Diam(indxDiamLessBreak(end):n)),1)];
    inx=[];
    iny=[];
    ix=newRBA.Diam(indxDiamLessBreak(end)).^(1:-1:0);
    iy=f_2s2fAll(i,indxDiamLessBreak(end));
    
    coefFit_2s2f_2All(i,:)=lsqlin(V,f_2s2fAll(i,indxDiamLessBreak(end):end),...
                              inx,iny,ix,iy);

    %evaluate coefficients obtained
    f_2s2fAll(i,indxDiamLessBreak(end):n)=coefFit_2s2f_2All(i,1)*...
                            newRBA.Diam(indxDiamLessBreak(end):n)+...
                            coefFit_2s2f_2All(i,2);
    res_2s2f=(f_2s2fAll(i,:)-newRBA.RBA')*(f_2s2fAll(i,:)-newRBA.RBA')';
    L_2s2f=-n/2*log(res_2s2f);
    BIC_2s2f(i)=L_2s2f-1/2*p_2s2f*log(n);
    
end

[~,indxmaxBIC22f]=max(BIC_2s2f);
 
breakpoint_2s2f=breakp(indxmaxBIC22f);


%% plot BIC of models

figure(1)
plot(breakp,BIC_1s*ones(numel(breakp),1),'LineWidth',2)
hold on
plot(breakp,BIC_2s,'LineWidth',2)
hold on
plot(breakp,BIC_2s2f,'LineWidth',2)
xlabel('Break point (km)')
ylabel('BIC')

legend('1 slope model','2 slopes model (1 fixed)',...
       '2 slopes model (2 free)','Location','northoutside')


%% plot data and models

figure(2)
% 1 slope
subplot(3,1,1)
plot(newRBA.Diam,newRBA.RBA,'.')
hold on
plot(newRBA.Diam,f_1s,'k','LineWidth',3)
ylabel('RBA (mGal)')
l1=legend('Data','Best fit line');
l1.Location='northeastoutside';
title('1 slope model with highest BIC')

%2 slopes 1 fixed
subplot(3,1,2)
plot(newRBA.Diam,newRBA.RBA,'.')
hold on
plot(newRBA.Diam,f_2sAll(indxmaxBIC2,:),'Color','k','LineWidth',3)
ylabel('RBA (mGal)')
hold on
xline(breakpoint_2s,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
l2=legend('Data','Best fit line','Breakpoint');
l2.Location='northeastoutside';
title('2 slopes model with highest BIC (1 slope fixed)')

%2 slopes 2 free

subplot(3,1,3)
plot(newRBA.Diam,newRBA.RBA,'.')
ylim([-100 100])
hold on
plot(newRBA.Diam,f_2s2fAll(indxmaxBIC22f,:),...
    'Color','k',...
    'LineWidth',3)
xlabel('Crater diameter (km)')
ylabel('RBA (mGal)')
hold on
xline(breakpoint_2s2f,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
l3=legend('Data','Best fit line','Breakpoint');
l3.Location='northeastoutside';
title('2 slopes model with highest BIC (2 slopes free)')

end