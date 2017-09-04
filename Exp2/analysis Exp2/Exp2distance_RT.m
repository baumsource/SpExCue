function [fig,result] = Exp2distance_RT(correct)
% evaluate response times of Exp2
% [fig,result] = Exp2distance_RT(correct)
% correct:  1 ... only correct trials (default)
%           0 ... all trials

if ~exist('correct','var')
  correct = 1;
end

datapath = fullfile('..');
expfn = {'distance'};
% expfn = expfn(session);

  expLabel = {''};

% tmp = load(fullfile(datapath,'analysis','SpExCue_Exp1eeg_subjects_tab.mat'));
% subject = tmp.subject.name;
% gender = tmp.subject.gender;
% sourceDirection = tmp.subject.direction;

subject = {'S07'};

conditions = {[0,1];[0,0.5];[0.5,1];[1,0];[0.5,0];[1,0.5]};
condLabelData = { '0-1';'0-.5';'.5-1';'1-0';'.5-0';'1-.5'};
condLabelPlot = { '0\leftrightarrow1';'0\leftrightarrow.5';'.5\leftrightarrow1';...
              };
response = {'receding','approaching'};

%% Percentages
pResp = nan(length(subject),length(expLabel),length(conditions),length(response));
seResp = pResp;
for ee = 1:length(expLabel)
  for ss = 1:length(subject)
    tmp = load(fullfile(datapath,'data',['Exp2behav_',subject{ss},'_',expfn{ee}]));
    E = tmp.subj.E;
    RT = tmp.subj.RT;
    Mcomb = tmp.subj.Mcomb;
    for cc = 1:length(conditions)
      idc = Mcomb(:,1) == conditions{cc}(1) & Mcomb(:,2) == conditions{cc}(2);
      N = sum(idc);
      RTidc = RT(idc);
      pResp(ss,ee,cc,1) = nanmedian(RTidc(E(idc) >  0)); % farther
      pResp(ss,ee,cc,2) = nanmedian(RTidc(E(idc) <  0)); % closer
      seResp(ss,ee,cc,1) = nanstd(RTidc(E(idc) >  0)); % farther
      seResp(ss,ee,cc,2) = nanstd(RTidc(E(idc) <  0)); % closer
    end
  end
end

% In percent
% pResp = 100*pResp; 

% Average constant conditions
% pResp = cat(3,pResp(:,:,1:6,:),mean(pResp(:,:,7:9,:),3));

% Quartiles
quartResp = quantile(pResp,[.25,.5,.75]); 

% Standard errors
% seResp = std(pResp,[],1)/sqrt(length(subject));

%% Percent correct

% % ANOVA
% RT = squeeze(cat(3,pResp(:,ee,1:3,1),pResp(:,ee,4:6,2))); % percent correct
% DV = array2table(RT); 
% contrast = strrep(condLabelPlot(1:3),'\leftrightarrow','-');
% contrast = repmat(contrast,[2,1]);
% direction = cell(length(contrast),1);
% direction(1:3) = {'approaching'};
% direction(4:6) = {'receding'};
% IVs = table(contrast,direction); %,'VariableNames',{'separation','direction'}
% rm = fitrm(DV,['RT1-RT',num2str(length(contrast)),' ~ 1'],'WithinDesign',IVs);
% [ranovaResult,~,C,~] = ranova(rm,'WithinModel','direction*contrast');
% 
% ranovaResult.Properties.RowNames = strrep(ranovaResult.Properties.RowNames,'(Intercept):','');
%     
% % Sphericity corrections
% spherCorr = epsilon(rm,C);
% 
% % Add corrected DFs to ranova table
% idrep = round(0.5:0.5:length(spherCorr.GreenhouseGeisser)); % repeat iteratively
% ranovaResult.DFGG = ranovaResult.DF .* ...
%   reshape(spherCorr.GreenhouseGeisser(idrep),size(ranovaResult.DF));
% 
% % Add effect sizes to ranova table
% SSeffect = ranovaResult.SumSq(1:2:end);
% SSerror = ranovaResult.SumSq(2:2:end);
% eta_pSq = nan(2*length(SSerror),1);
% eta_pSq(1:2:end) = SSeffect./(SSeffect+SSerror); % effect size per (eta_partial)^2
% ranovaResult.eta_pSq = eta_pSq;
% 
% disp(ranovaResult)
% 
% multcompare(rm,'contrast')


%% Plot
fig = figure;
ha = tight_subplot(1,length(expLabel),0,[.13 .05],[.13 .01]);
x = [1:3,1:3];
XTick = x(1:3);
XLim = [0.4,3.6];
dx2 = .05*[-1,0,1]; % between responses
symb = '^v';
lineStyle = {':';'-';''}; % increase, decrease
color = [41,17,97;149,123,109;229,68,0]/255;
color = color([1,3],:);

if correct
  pResp(:,:,4:6,1) = nan;
  pResp(:,:,1:3,2) = nan;
end

for ee = 1:length(expLabel)
  axes(ha(ee))

  for rr = 1:length(response)
    y = nanmean(pResp(:,ee,:,rr),1);
    l = seResp(1,ee,:,rr);
    u = seResp(1,ee,:,rr);
    idx = {1:3;4:6};
    for ii = 1:length(idx)
      h(rr,ii) = errorbar(x(idx{ii})+dx2(rr),y(idx{ii}),l(idx{ii}),u(idx{ii}),...
        [symb(rr),lineStyle{ii}]); 
      hold on
    end
    set(h(rr,:),'MarkerFaceColor',color(rr,:),'Color',color(rr,:))
  end
  set(h(:,1),'MarkerFaceColor','w')
  set(gca,'XTick',XTick,'XTickLabel',condLabelPlot)

  title(expLabel{ee})
  if ee==1
    ylabel('Response Time (s)')
    XLabel = 'Spectral contrast switch, C_{1,2}';
%     if length(session) == 2
%       XLabel = [repmat(' ',1,50),XLabel];
%     end
    xlabel(XLabel)
  else
    set(gca,'YTickLabel',{})
  end
end

%% Plot legend
if correct
  legend(h([1,4]),response)
else

  fig(2) = figure;
  direction = {'C_1<C_2','C_1\geqC_2'};
  yspacing = 0.75;
  ymax = yspacing*(length(response)+1);
  dxResp = 1.2;
  for ii = 1:2
    for rr = 1:length(response)
      y = yspacing*(4-rr);
      if ii == 1
        t(ii,rr) = text(ii-dxResp,y,response{rr},'HorizontalAlignment','center');
      end
      hold on
      l(ii,rr) = plot(ii+[-.4,.4],[y,y],lineStyle{ii});
      sy(ii,rr) = plot(ii,y,[symb(rr),lineStyle{ii}]);
      set([l(ii,rr),sy(ii,rr)],'Color',color(rr,:),'MarkerFaceColor',color(rr,:))
    end
    text(ii,ymax,direction{ii},'HorizontalAlignment','center')
  end
  text(0,ymax,'response','HorizontalAlignment','center')
  set(sy(1,:),'MarkerFaceColor','w')
  set(gca,'XTick',[],'YTick',[])
  axis([1-dxResp-yspacing,2+yspacing,1-yspacing,ymax+0.9*yspacing])
  box off
  % axis equal

  copyobj(get(fig(2),'Children'),fig(1))
  ax = get(fig(1),'Children');
%   if length(session) == 2
%     set(ax(1),'Position',[0.45    0.45    0.43    0.17])
%   else
    set(ax(1),'Position',[0.54    0.17    0.45    0.30]) % y=0.67 for top alignment
    set(ax(2),'Position',[0.12    0.17    0.42    0.80])
%   end

end

%% Output
result = [];
for ee = 1:length(expLabel)
  result.exp(ee).name = expLabel{ee};
  for ss = 1:length(subject)
    result.exp(ee).subj(ss).name = subject{ss};
%     result.exp(ee).subj(ss).gender = gender{ss};
%     result.exp(ee).subj(ss).direction = sourceDirection{ss};
%     result.exp(ee).subj(ss).bias = bias(ss,ee);
    result.exp(ee).subj(ss).RT = array2table(squeeze(pResp(ss,ee,:,:)),...
      'RowNames',condLabelData,'VariableNames',response);
  end
end
end