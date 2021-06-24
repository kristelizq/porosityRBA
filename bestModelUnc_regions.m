
function [bestCoef1High,bestCoef2High,bestBreakHigh,...
          upperBHigh,lowerBHigh,...
          upperS1High,lowerS1High,...
          upperS2High,lowerS2High,...
          bestCoef1SPA,bestCoef2SPA,bestBreakSPA,...
          upperBSPA,lowerBSPA,...
          upperS1SPA, lowerS1SPA,...
          upperS2SPA, lowerS2SPA,...
          bestCoeffMare,...
          upperSmare, lowerSmare]=...
          bestModelUnc_regions(nsamples,printing)
%This function finds the best fitting model and uncertainty 
%for RBA data of highlands, spa and mare regions.
%For highlands and spa: 2 slope models with 2 free slopes.
%For mare: 1 slope model.
%This function takes ~1 hour to run

%% load RBA data and separate by region

%RBA data
newRBA=readtable('../../../data/RBArobbins18Goossens20_taper.csv');
craterData=readtable('../../../data/regionCraters.csv');

%sort RBA by diameter 
[~,indxSort]=sort(newRBA.Diam);
newRBA=newRBA(indxSort,:);
craterData=craterData(indxSort,:);

indxMare=find(craterData.mare);
indxSPA=find(craterData.insSPA);
indxHighlands=find(craterData.mare==0 & craterData.insSPA==0);

rbaMare=newRBA(indxMare,:);
rbaSPA=newRBA(indxSPA,:);
rbaHighlands=newRBA(indxHighlands,:);

%% calculate best fit models and uncertainty for each region

%highlands
[bestCoef1High,bestCoef2High,bestBreakHigh,bestFHigh,...
 upperBHigh,lowerBHigh,upperS1High,lowerS1High,...
 upperS2High,lowerS2High,...
 xCoorFHigh,yCoorFHigh]=model2s2fAux(rbaHighlands,nsamples);

%spa
[bestCoef1SPA,bestCoef2SPA,bestBreakSPA,bestFSPA,...
 upperBSPA,lowerBSPA,...
 upperS1SPA,lowerS1SPA,...
 upperS2SPA,lowerS2SPA,...
 xCoorFSPA,yCoorFSPA]=model2s2fAux(rbaSPA,nsamples);

%mare
[bestCoeffMare,bestFmare,upperSmare,...
 lowerSmare,xCoorFmare,yCoorFmare]=model1sAux(rbaMare,nsamples);

%% Plot best fit models and uncertainty of highlands, mare and spa

figure(1)

%highlands----------------------------------------------------
subplot(3,1,1)
plot(rbaHighlands.Diam,rbaHighlands.RBA,'.')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorFHigh,yCoorFHigh,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(rbaHighlands.Diam,bestFHigh,...
    'Color','k',...
    'LineWidth',3)
hold on
xline(bestBreakHigh,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
hold on
xline(lowerBHigh,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
hold on
xline(upperBHigh,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
legend('Data','95% CI slope','Best fit line',...
       'Breakpoint','95% CI breakpoint');

title('Highlands')
pbaspect([2,1,1])
yl = ylim;
ylim(yl)


%spa----------------------------------------------------
subplot(3,1,2)
plot(rbaSPA.Diam,rbaSPA.RBA,'.')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorFSPA,yCoorFSPA,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(rbaSPA.Diam,bestFSPA,...
    'Color','k',...
    'LineWidth',3)
hold on
xline(bestBreakSPA,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
hold on
xline(lowerBSPA,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
hold on
xline(upperBSPA,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
title('SPA')
pbaspect([2,1,1])
ylim(yl)

%mare-----------------------------
subplot(3,1,3)
plot(rbaMare.Diam,rbaMare.RBA,'.')
xlabel('Crater diameter (km)')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorFmare,yCoorFmare,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(rbaMare.Diam,bestFmare,...
    'Color','k',...
    'LineWidth',3)
title('Mare')
pbaspect([2,1,1])
ylim(yl)

% %move legend out of the way
splots=get(gcf,'children');
set(splots(3),'position',[0.7 0.72 0.2 0.169]);

%% print figure
%Print figure
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'white')
posfig=get(gca,'Position');
set(gca,'Position',posfig);
set(gcf,'Units','Inches');
pos = get(gcf,'position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if printing=='yes'
    filenamef=['modelAndUnc_regions' ];
    print(gcf,'-painters','-dpdf',filenamef);
end


%zoom in ----------------------------
%-----------------------
figure(2)

%highlands----------------------------------------------------
subplot(3,1,1)
plot(rbaHighlands.Diam,rbaHighlands.RBA,'.')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorFHigh,yCoorFHigh,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(rbaHighlands.Diam,bestFHigh,...
    'Color','k',...
    'LineWidth',3)
hold on
xline(bestBreakHigh,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
hold on
xline(lowerBHigh,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
hold on
xline(upperBHigh,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
% legend('Data','95% CI slope','Best fit line',...
%        'Breakpoint','95% CI breakpoint');

title('Zoom in Highlands')
pbaspect([2,1,1])
ylim([-10,10])


%spa----------------------------------------------------
subplot(3,1,2)
plot(rbaSPA.Diam,rbaSPA.RBA,'.')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorFSPA,yCoorFSPA,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(rbaSPA.Diam,bestFSPA,...
    'Color','k',...
    'LineWidth',3)
hold on
xline(bestBreakSPA,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
hold on
xline(lowerBSPA,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
hold on
xline(upperBSPA,...
     'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2)
title('Zoom in SPA')
pbaspect([2,1,1])
ylim([-10,10])

%mare-----------------------------
subplot(3,1,3)
plot(rbaMare.Diam,rbaMare.RBA,'.')
xlabel('Crater diameter (km)')
ylabel('RBA (mGal)')
hold on
p=patch(xCoorFmare,yCoorFmare,'g');
p.EdgeColor='g';
p.FaceAlpha=0.5;
p.EdgeAlpha=0.5;
hold on
plot(rbaMare.Diam,bestFmare,...
    'Color','k',...
    'LineWidth',3)
title('Zoom in Mare')
pbaspect([2,1,1])
ylim([-10,10])

%% print figure
%Print figure
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'white')
posfig=get(gca,'Position');
set(gca,'Position',posfig);
set(gcf,'Units','Inches');
pos = get(gcf,'position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if printing=='yes'
    filenamef=['modelAndUnc_regionsZoomIn' ];
    print(gcf,'-painters','-dpdf',filenamef);
end

end