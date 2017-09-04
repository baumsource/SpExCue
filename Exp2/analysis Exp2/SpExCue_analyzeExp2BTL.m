function SpExCue_analyzeExp2BTL(ID,fnext,printFlag)
% SpExCue_analyzeExp2BTL(ID,fnext,printFlag)

if not(exist('printFlag','var'))
    printFlag = false;
end
if not(exist('OptiPt','file'))
    addpath(fullfile('..','..','MATLAB_general'))
end

tmp = load(fullfile('..','data',['Exp2behav_' ID,'_',fnext]));
subj = tmp.subj;
E = subj.E; % closer
M = unique(subj.Mcomb);
NM = length(M);
iM = repmat(1:NM,[NM,1]);
iMtransp = transpose(iM);
iM = [iMtransp(:),iM(:)];
Mcomb = M(iM); % combinations of M
Ncomb = size(Mcomb,1);
%%
nResp = nan([NM,NM]);
for cc = 1:Ncomb
  idc = subj.Mcomb(:,1) == Mcomb(cc,1) & subj.Mcomb(:,2) == Mcomb(cc,2);
%   N = sum(idc);
  nResp(cc) = sum(E(idc)<0); % closer
end

%% Run BTL model fit
% nResp needs to be ordered such that first column of Mcomb is increasing
% first. Hence, row is M1 col is M2
disp('RB: BTL model fitted on "closer" responses')
A = {1;2;3;4}; % indices
Alabel = {'C = 0','C = .5','C = 1','KEMAR'};
[Ebtl,chistat] = OptiPt(nResp,A);
disp(chistat)

Ebtl = Ebtl/max(Ebtl); % normalize

%% Sort (important for KEMAR)
[Esort,isort] = sort(Ebtl);

%% Plot
fig = figure; 
plot(Esort,'.-')
set(gca,'XTick',1:NM,'XTickLabel',Alabel(isort),'XLim',[0.5,NM+0.5])
ylabel('Externalization')
title([ID,' at ',num2str(subj.pos(1,1)),' deg'])

%% Save
if printFlag
  FontSize = 8;
  Resolution = '-r600';
  LineWidth = 1;
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(findall(fig,'-property','LineWidth'),'LineWidth',LineWidth)
  FigSize = [8,8];
  set(fig,'PaperUnits','centimeters','PaperSize',FigSize,'PaperPosition',[0,0,FigSize])
  print(fig,Resolution,'-dpng',fullfile(mfilename,[mfilename,'_',ID]))
end