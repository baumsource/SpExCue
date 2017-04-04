function h = SpExCue_corrPlot(x,y)
% SpExCue_corrPlot - scatter plot with regression line and R-square value


% Evaluate regression
tbl = table(x(:),y(:),'VariableNames',{'x','y'});
lme = fitlme(tbl,'y~x');

% Axis limits and text coordinates + correlation sign
dr = 0.1;
dx = dr*range(x);
dy = dr*range(y);
Lim = [min(x)-dx,max(x)+dx,min(y)-dy,max(y)+dy];
if lme.Coefficients.Estimate(2) > 0 % R2 left 
  xtxt = Lim(1)+dx;
  rsign = '';
else
  xtxt = Lim(2)-4*dx;
  rsign = '-';
end

% p-value stars
if lme.coefTest < .001
  pvaltxt = '***';
elseif lme.coefTest < .01
  pvaltxt = '**';
elseif lme.coefTest < .05
  pvaltxt = '*';
else
  pvaltxt = '';
  disp(lme.coefTest)
end

rtxt = ['r = ',rsign,num2str(sqrt(lme.Rsquared.Ordinary),'%1.2f'),pvaltxt];
  
% Plot
h(2) = plot(x,lme.fitted,'k-');
hold on
h(1) = plot(x,y,'ok');
set(h(1),'MarkerFaceColor','w')
text(xtxt,Lim(4)-dy,rtxt)

axis(Lim)

end