function [h,p,stats] = SB_permutationTest(a,b,N,thresh,alpha,tail)
%
% [h,p,stats] = SB_permutationTest(a,b,N,thresh,alpha,tail)
%   Performs cluster-based nonparametric statistics on two times series
%   across two different conditions.  Assumes same subjects in both
%   conditions.  The null distribution of the test statistic is derived
%   from N bootstrapped permutations of the data.
%
%   The test statistic used in this version is the sum of the t-statistic
%   within a time cluster that exceeds a certain thresholds, thresh,
%   defined by the user.  The choice of the threshold value is arbitrary,
%   but does have implications in terms of the sensitivity of the test.
%   It does not affect False Alarm rate of the test.
%
%   For a more detailed description of the methods and statistical validity
%   see: 
%
%   Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing 
%       of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 
%       177?190. http://doi.org/10.1016/j.jneumeth.2007.03.024
%
% INPUT VARIABLES
%       a : multi-subject time series for Condition A <subj x samples>
%       b : multi-subject time series for Condition B <subj x samples>
%       N : number of bootstrap permutations (recommend at least 1000)
%  thresh : user selected t-statistic threshold (default=0.95)
%   alpha : p-value criteria
%    tail : direction of the one-sided t-test ['left' or 'right']
%
% OUTPUT VARIABLE
%       h : sample-by-sample test decision variable 
%               1=reject null hypothesis at significance level, alpha
%               0=failure to reject null hypothesis
%       p : cluster-by-cluster p-values
%   stats : structure of test results
%
% Created by: Scott Bressler 08-Jun-2016 (last modified by Robert 12-Jul-2016)

if nargin < 4
    thresh = 0.95; % t-distribution threshold
    alpha = 0.05;  % significance level
    tail = 'left'; % t-test tail
elseif nargin < 5
    alpha = 0.05;
    tail = 'left';
elseif nargin < 6
    tail = 'left';
end

% Initialize test decision variable, h
h = zeros(1,size(a,2));

% Determine number of subjects per condition
nA = size(a,1);
nB = size(b,1);
df = nA-1;  % degrees of freedom

% % Determine correct t-stat threshold based on tail
% if(strcmp(tail,'right'))
%     THR = tinv(thresh,df);
% elseif(strcmp(tail,'left'))
%     THR = tinv(1-thresh,df);
% end

% Id = cat(1,ones(nA,1),2*ones(nB,1)); % delete later

% Calculate the "Observed" test statistic
[Xobs.h,Xobs.p,Xobs.ci,Xobs.stats] = ttest(a,b,'Tail',tail);
[ObsClstrStat,kClstr] = findClusters(Xobs.stats.tstat,thresh,df,tail);

% Return null results if no clusters in Observed Data
if(isempty(ObsClstrStat));
    h = zeros(1,size(a,2));
    Ci = [];
    p = [];
    stats = struct('ClusterStats',[],'ClusterIndices',Ci,'pValue',[],...
                   'p',p,'threshold',thresh,'alpha',alpha,...
                   'numPermutations',N,'method',sprintf('t-Tailed:%s',tail));
    fprintf('WARNING: No clusters in OBSERVED DATA.\n');
    fprintf('         Consider adjusting threshold variable\n');
    return;
end

%% Prepare data for bootstrap permutation procedure
Xcat = cat(1,a,b); % gather both groups into a single matrix

prmStat = zeros(N,1); % initialize variable

for k = 1:N % loop through all bootstrap trials
%     Ik = randperm(nA+nB); % randomly permute subject data
    Ia = rand(nA,1) > 0.5;
    Ia = cat(1,Ia,not(Ia));
    
    % Calculate permutation test statistic (difference a-b)
%     [~,~,~,Xtest] = ttest(Xcat(Ik(1:nA),:),Xcat(Ik(nA+1:end),:),'Tail',tail);
    [~,~,~,Xtest] = ttest(Xcat(Ia,:),Xcat(not(Ia),:),'Tail',tail);
    
    % Locate time serires clusters
    [clstr,~] = findClusters(Xtest.tstat,thresh,df,tail);
    
    % Determine max/min cluster-level permuted test statistic
    if ~isempty(clstr)
        if(strcmp(tail,'left'))
            prmStat(k) = min(clstr);
        elseif(strcmp(tail,'right'));
            prmStat(k) = max(clstr);
        end
    end
    
end

% Calculate permuation null distribution, based on empirical cumulative
% distribution function
[f,x] = ecdf(prmStat);

% Determine alpha significance value, pVal, and test decision variable
switch tail
    case 'right'
        pVal = x(find(f>=(1-alpha),1,'first'));
        if any(ObsClstrStat>=pVal)
            for n = find(ObsClstrStat>=pVal)
                h(kClstr{n}) = 1;
            end
        end
    case 'left'
        pVal = x(find(f<=alpha,1,'last'));
        if any(ObsClstrStat<=pVal)
            for n = find(ObsClstrStat<=pVal)
                h(kClstr{n}) = 1;
            end
        end
end


% Determine cluster-level p-values
p = ones(length(ObsClstrStat),1);
for n = 1:length(ObsClstrStat)
    switch tail
        case 'right'
            if ObsClstrStat(n)>max(x)
                p(n) = 0;
            else
                p(n) = 1-f(find(x>=ObsClstrStat(n),1,'first'));
            end
        case 'left'
            if ObsClstrStat(n)<min(x)
                p(n) = 0;
            else
                p(n) = f(find(x>=ObsClstrStat(n),1,'first'));
            end
    end
end

% Output stats summary structure
stats.ClusterStats = ObsClstrStat;
stats.ClusterIndices = kClstr;
stats.pValue = pVal;
stats.p = p;
stats.threshold = thresh;
stats.alpha = alpha;
stats.numPermuations = N;
stats.method = sprintf('t-Tailed:%s',tail);
        
end

function [clStat,CL] = findClusters(Xstat,thresh,df,tail)

if(strcmp(tail,'right'))
    cl = Xstat >= tinv(thresh,df);
elseif(strcmp(tail,'left'))
    cl = Xstat <= tinv(1-thresh,df);
elseif(strcmp(tail,'both'))
    cl = abs(Xstat) >= tinv(thresh/2,df);
end
        
cl = [0;cl(:)]'; % zero pad for initial condition

clStat = [];
CL = [];
clNum = 0;
temp = [];

for k = 2:length(cl)
    if(cl(k))
        temp = cat(2,temp,k);
    end
    
    if(cl(k)==0 && cl(k-1)==1)
        clNum = clNum+1;
        CL{clNum} = temp-1;
        temp = [];
    end
end

if ~isempty(CL)
    for k = 1:length(CL)
        clStat(k) = sum(Xstat(CL{k}));
    end
end


% if cl(end) % Special case where last sample is a single point "cluster"
%     cl(end) = 0; % get rid of it
% end
% 
% cl1 = find(cl==1);
% cl0 = find(cl==-1)-1;
% 
% if isempty(cl1)
%         clStat = [];
%         CL = [];
%         return;
% end
% 
% %% Address special clustering cases
% if(cl0(1)<cl1(1))
%     cl0 = cl0(2:end);
% end
% 
% if(cl0(end)<cl1(end))
%     cl0(end) = length(cl);
% end
% 
% %%
% 
% 
% CL = [cl1',cl0'];
% for n = 1:size(CL,1)
%     clStat(n) = sum(Xstat(CL(n,1):CL(n,2)));
% end

end
    