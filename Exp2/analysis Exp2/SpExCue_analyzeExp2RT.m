function SpExCue_analyzeExp2RT(ID,fnext,printFlag)
% SpExCue_analyzeExp2RT - Analysis of behavioral data 
%  Usage: SpExCue_analyzeExp2RT(ID,fnext,printFlag)
%  Input parameters:
%  ID:        subject ID 
%  fnext:     experimental postfix, e.g., LR or distance
%  printFlag: set true to plot results (default: true)
%
%  Output parameters: 

if not(exist('printFlag','var'))
    printFlag = false;
end
if not(exist('OptiPt','file'))
    addpath(fullfile('..','..','MATLAB_general'))
end

%% Load data
tmp = load(fullfile('..','data',['Exp2behav_' ID,'_',fnext]));
subj = tmp.subj;
RT = 1000*subj.RT; % RT in ms
M = unique(subj.M);
NM = length(M);
if NM ==3
    iorder = [1,3,2];
    lbl = {'flat','KEMAR','individual'};
    M = M(iorder);
end
lat = unique(subj.pos(:,1));
%%
mRT = nan(NM,length(lat));
q1RT = nan(NM,length(lat));
q2RT = nan(NM,length(lat));
for mm = 1:NM
    for ll = 1:length(lat)
        idc = subj.M == M(mm) & subj.pos(:,1) == lat(ll);
        disp(sum(idc))
        mRT(mm,ll) = median(RT(idc));
        q1RT(mm,ll) = prctile(RT(idc),25);
        q2RT(mm,ll) = prctile(RT(idc),75);
    end
end

%% Plot
fig = figure; 
XTick = 1:NM;
dx=0.05;
errorbar(XTick-dx,mRT(:,1),q1RT(:,1),q2RT(:,1),'ro-')
hold on
errorbar(XTick+dx,mRT(:,2),q1RT(:,2),q2RT(:,2),'bo-')
set(gca,'XTick',XTick,'XTickLabel',lbl,'XLim',[0.5,NM+0.5])
ylabel('Response time (ms)')
title([ID,' at ',num2str(subj.pos(1,1)),' deg'])
legend('right','left')

%% Save
if printFlag
  FontSize = 8;
  Resolution = '-r600';
  LineWidth = 1;
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(findall(fig,'-property','LineWidth'),'LineWidth',LineWidth)
  FigSize = [8,8];
  set(fig,'PaperUnits','centimeters','PaperSize',FigSize,'PaperPosition',[0,0,FigSize])
  fn = fullfile(mfilename,[mfilename,'_',ID,'_',fnext]);
  print(fig,Resolution,'-dpng',fn)
end