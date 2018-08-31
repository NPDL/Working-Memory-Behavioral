% WMSpan Project - Bedny Lab %%%%%%%%%%%%%%%
% 2018
% Data collected by Rita Loiotile
% Data analyzed by Rita Loiotile, Karen Arcos, and Nora Harhen
% Script by Nora Harhen (nharhen1@jhu.edu) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze complex span task data on the main task and on the distractor
% task

clear all;
load('allData_28-Aug-2018 11:27:13.mat')
%% Load and combine data 
subjects{1} = [1:10 12:16 18:23]; %S subjects, eliminated o w/o critical tasks
subjects{2} = [1:2 4:5 7:10 13:16 18:25]; %CB subjects, eliminated ones w/o critical tasks
N = cellfun(@(x) length(x), subjects); % find number of subjects in each group 
%columns we care about for this analysis
toKeep = [13 15]; % proportion correct positions, prop corr letters
toKeepEqn = [8,10,14:15];
headingsToKeep = headings{1}(toKeep); %headings are same for both tasks
%Our two groups 
vision = {'S', 'CB'}; 
span = [2:9]; % highest span someone went to was 9
% Intialize matrices
mathPerfperSpan = zeros(2,25,9)+0.5; % First dim: group, second dim: subject, third dim: span
mathPerfAc = zeros(2,25,2)+0.5;% First dim: group, second dim: subject, third dim: incorrect or correct
mathPerfSpanDepen = zeros(2,25,2)+0.5;% First dim: group, second dim: subject, third dim: incorrect or correct
%create matrix of data
for v = 1:length(subjects) % loop through our two groups
    for su = subjects{v} %loop through subjects in groups
        clear tempMat
        %create temporary matrix
        tempMat = eval(['cell2mat(',vision{v},'Data(',num2str(3),').su',num2str(su),'.recall(:,toKeep))']);
        if size(tempMat,1) < 24 % need to zero pad if reached less than span 9
            numRowsAdd = 24 - size(tempMat,1);
            rowsToAdd= [0 0];
            tempMat = [tempMat;repmat(rowsToAdd,numRowsAdd,1)];
        end
        %average across span, each span has two trials associated with
        %it, so we can average every other row
        n=3;
        averages = arrayfun(@(i) mean(tempMat(i:i+n-1,:)),... % averages is a cell array
            1:n:length(tempMat)-n+1,'UniformOutput', false);
        averages = cell2mat(averages); % convert from cell to matrix
        %reshape data from vector to Nx2x8 first dim: subject, second dim:
        %measure, third dim: span
        %save to giant data matrix
        eval(['intraSubAvg.', vision{v},'(',num2str(su),',:,:) = reshape(averages,[2,8]);']);
        
        % get the % correct performance on the math equations per span
        cellMat = eval([vision{v},'Data(',num2str(3),').su',num2str(su),'.eqtn(:,toKeepEqn);']);
        mathPerf = cellfun(@(x) double(x),cellMat);
        % need to convert the third column from having trialNum to having
        % spanInfo 
        startSpanInd = 1:3:mathPerf(end,3);
        for i =1:length(startSpanInd)
            ind = ismember(mathPerf(:,3),startSpanInd(i):startSpanInd(i)+2);
            mathPerf(ind,3) = span(i);
        end
        % average over trials per 
        spanEqn = [0 2:max(mathPerf(:,3))]; 
        for sp = 1:numel(spanEqn)
            spanOfInterest = mathPerf(:,3) == spanEqn(sp);
            mathPerfperSpan(v,su,sp) = mean(mathPerf(spanOfInterest,2));
        end
        % get the performance on math if correct or incorrect on that main
        % task trial  and if the eqtn itself is correct or not 
        for a = [0 1]
            accOfInterest = mathPerf(:,1) == a; spanAccOfInterest = mathPerf(:,end) == a;
            mathPerfAc(v,su,a+1) = mean(mathPerf(accOfInterest,2)); 
            mathPerfSpanDepen(v,su,a+1) = mean(mathPerf(spanAccOfInterest,2)); 
        end
    end
end 


%% Set default settings for plotting 
% Because i am lazy and don't wanna set this every time
set(groot,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTitleFontSize',1)
set(groot,'DefaultLineMarkerSize',10)
set(groot,'DefaultLineLineWidth',2)
color = [244 140 66;66 119 244]./255; %s is first row (orange) , cb is second row (blue)

%% Performance On Main Task 

%%% Across Spans %%% 
[statPerSpan(1)] = descStats(intraSubAvg.S(subjects{1},:,:));
[statPerSpan(2)] = descStats(intraSubAvg.CB(subjects{2},:,:));

% Plot the data %
figure;
for plt = 1:2 % loop through the two measures (letters v positions)
    subplot(1,2,plt)
    hold on
    for v = 1:2 % loop through the two groups    
        points = squeeze(statPerSpan(v).mean(:,plt,:));
        errBars = squeeze(statPerSpan(v).ci(:,plt,:));
        
        e(v) = errorbar(points,errBars, 'o'); % plot the marker and errorbar
        e(v).Color = color(v,:); % set color of errorbars
        e(v).MarkerFaceColor = color(v,:);  % set color of markers
        e(v).MarkerSize = 10;  % set size of markers
    end % for v
    set(gca,'xtick',1:8,'xticklabels',span) % list the spans
    ylabel('% Correct'); xlabel('Span'); title(headingsToKeep{plt})
    ylim([0 1])
    legend('S','CB')
end

% Statistical Testing %
[h,p,~,stats] = ttest2(intraSubAvg.S(subjects{1},:,:),intraSubAvg.CB(subjects{2},:,:));

%%% Overall Performance %%%

% Calc Stats % 
[statOverall(1)] = descStats(mean(intraSubAvg.S(subjects{1},:,:),3));
[statOverall(2)] = descStats(mean(intraSubAvg.CB(subjects{2},:,:),3));

% Plot Data % 
measures = {'Positions','Letters'};
figure;
hold on;
for v = 1:2 %loop through group
    % select data to plot
    marker = squeeze(statOverall(v).mean);
    errBar = squeeze(statOverall(v).ci);
    e(v)= errorbar(marker,errBar,'o-', 'MarkerFaceColor',color(v,:));
    % plot data
    [e(v).Color, e(v).MarkerFaceColor] = deal(color(v,:));
end
ylabel('Proportion correct');
xlim([0.5 2.5]); ylim([0.25 0.75])
set(gca,'xtick',1:2,'xticklabels',measures)
set(gca, 'fontsize',16)
legend('S', 'CB')
title('Overall Performance')

% Statistical Testing % 
[h,p,~,stats] = ttest2(mean(intraSubAvg.S(subjects{1},:,:),3),mean(intraSubAvg.CB(subjects{2},:,:),3));
disp('Did t-test between groups on overall performance. Found signifigant differences on :')
disp(measures(find(h==1)))
disp(['P-value: ', num2str(p(h==1)), ' T-value: ', num2str(stats.tstat(h==1))])

%% Performance on Equations 

%%% Calculate Stats %%%
measure = {'PerfperSpan','PerfAc','PerfSpanDepen'};
for m = 1:numel(measure)
    for v = 1:2
        eval(['[stat',measure{m},'(v)] = descStats(squeeze(math',measure{m},'(v,subjects{v},:)));']);
    end
end

%%% Plot Results %%%

% Performance separated by Span
figure
hold on;
for v = 1:2
    % grab the data to plot
    markers = squeeze(statPerfperSpan(v).mean(:,:,:));
    errBars = squeeze(statPerfperSpan(v).se(:,:,:));
    %plot the data
    e(v) = errorbar(markers,errBars, 'o'); % plot the marker and errorbar
    [e(v).Color,e(v).MarkerFaceColor] = deal(color(v,:)); % set color of markers and errorbars
end
set(gca,'xtick',1:9,'xticklabels',[0 span]) % list the spans
ylabel('% Correct'); xlabel('Span'); title('Performance on Eqtns per Span')
ylim([0.5 1])
legend('S','CB')

% Performance separated by correct or incorrect equation overlaid with
% Performance separated by whether that main task trial was correct or
% incorrect 
lineStyles = {'o-', 'o--'}; % solid line for positions and dashed for letters 
lineColors = {'r','b'}; % red line for sighted and blue line for blind

figure
hold on;
for m = 2:3
    for v = 1:2
        % grab the data to plot
        strToeval = ['stat',measure{m},'(v).mean(:,:,:)']; 
        markers = squeeze(eval(strToeval));
        strToeval = ['stat',measure{m},'(v).se(:,:,:)']; 
        errBars = squeeze(eval(strToeval));
        %plot the data
        e(v) = errorbar(markers,errBars, [lineColors{m-1}, lineStyles{m-1}], 'color' ,color(v,:)); % plot the marker and errorbar
        %[e(v).Color,e(v).MarkerFaceColor] = deal(color(v,:)); % set color of markers and errorbars
    end
end
set(gca,'xtick',1:2,'xticklabels',{'Correct', 'Incorrect'}) % list the spans
ylabel('% Correct'); xlabel('Type'); 
ylim([0.7 1]); xlim([0.5 2.5])
legend('S Eqtn','CB Eqtn', 'S Trial', 'CB Trial')


%%% Statistical Testing %%%
allSpans = [0 span];
[h,p,~,stats] = ttest2(mathPerfperSpan(1,subjects{1},:), mathPerfperSpan(2,subjects{2},:));
disp(['Did a between group t-test on equation per span. These spans are signfigant:',num2str(allSpans(find(h==1)))])
disp(['P-value: ' num2str(p(find(h==1))), ' T-value: ', num2str(stats.tstat(find(h==1)))])
% Performance separated by whether the equation was correct or incorrect
conditions = {'Correct Eqtn', 'Incorrect Eqtn'};
[h,p,~,stats] = ttest2(mathPerfAc(1,subjects{1},:), mathPerfAc(2,subjects{2},:));
disp('Did a between group t-test on equation depending on whether eqtn was correct or incorrect. ')
disp(['The ', conditions{find(h==1)},' condition was signfigant:'])
disp(['P-value: ' num2str(p(find(h==1))), ' T-value: ', num2str(stats.tstat(find(h==1)))])
% Peformance depending on whether or not they got that trial correct or
% incorrect 
conditions = {'Correct Trial', 'Incorrect Trial'};
[h,p,~,stats] = ttest2(mathPerfSpanDepen(1,subjects{1},:), mathPerfSpanDepen(2,subjects{2},:));
disp('Did a between group t-test on equation depending on whether eqtn was correct or incorrect. ')
disp(['The ', conditions{find(h==1)},' condition was signfigant:'])
disp(['P-value: ' num2str(p(find(h==1))), ' T-value: ', num2str(stats.tstat(find(h==1)))])

%% Linear Regression
% predict performance with predictors: sight, span, accuracy on equations,
% sight*span, sight*accuracy equations, span*accuracyeqtns, triple
% interatcion

%%% Lump Groups Together %%%

linearRegressNames = {'Intercept' ,'isBlind', 'Span', 'Accuracy on Equations',...
    'isBlind*Span','isBlind*AccEqtns','Span*AccEqtns','isBlind*Span*AccEqtns'};
% let's start with proportion of positions correct 
% intialize our y hat vector and our regressor vector
Y = nan(sum(N)*length(span),1);
X = nan(sum(N)*length(span),7);
allSubs = [subjects{1} subjects{2}];
startSubPerf = 1:length(span):sum(N)*length(span); % where to start filling for each subject
% Fill our y hat and regressor vectors 
for su = 1:length(allSubs)
        blind = (su > N(1)) + 1; % 1 = s, 2 = cb
        groupSubNum = allSubs(su); 
        subjectRange = startSubPerf(su):startSubPerf(su)+7; 
        for load = span
            ind = subjectRange(load-1); 
            Y(ind) = eval(['squeeze(intraSubAvg.' vision{blind},'(groupSubNum,1,load-1))']); % we want to predict average performance
            X(ind,1:3) = [blind-1,load,mathPerfperSpan(blind,groupSubNum,load)];
            X(ind,4:7) = [X(ind,1)*X(ind,2), X(ind,1)*X(ind,3),X(ind,2)*X(ind,3), X(ind,1)*X(ind,2)*X(ind,3)];
        end
end

% z-score predictors so in standard units
X = zscore(X); 
% ourglm
[bLin,dev,stats] = glmfit(X,Y);

%%% Separate by Groups %%%

% Re-do linear regression by separate by group and get beta weights for each
% subject
linGroupRegressorNames = {'Intercept','Span', 'Accuracy on Equations', 'Span*AccEqtns'};
for v = 1:2
    for su = subjects{v}
        clear X Y
        Y = nan(8,1);
        X = nan(8,3);
        for load = 1:length(span)
            Y(load) = eval(['squeeze(intraSubAvg.', vision{v},'(su,1,load));']);
            X(load,1:2) = [load mathPerfperSpan(v,su,load+1)'];
        end
        % add interaction
        X(:,3) = X(:,1).*X(:,2);
        % z-score X
        X = zscore(X);
        % do regression and save results for each group
        [bLinGroup(v,su,:),~,stats(v,su,:)] = glmfit(X,Y);
    end
end

%%% Statistical Testing %%%

% Test each group against zero
for v=1:2
    clear h p
    [h,p] = ttest(bLinGroup(v,subjects{v},:));
    disp(['These regressors are signifigant for the ', vision{v}, ' group'])
    disp(linGroupRegressorNames(find(h==1)))
    disp(['P-values: ', num2str(p(h==1))])
end

% Do T-Test between groups
clear h p 
[h,p] = ttest2(bLinGroup(1,subjects{1},:),bLinGroup(2,subjects{2},:)); 
disp('Beta weights signifigantly differ between the groups for these regressors:')
disp(linGroupRegressorNames(find(h==1)))


% check for correlation between span performance and eqtn at same span
% performance
for v = 1:2
    X = squeeze(eval(['intraSubAvg.',vision{v},'(subjects{v},1,:)']));
    Y = squeeze(mathPerfperSpan(v,subjects{v},2:9));
    [rho,p] = corr(X,Y,'type','spearman'); % we only care about the correlations between the same span
    rhoDiag(v,:) = diag(rho)'; pDiag(v,:) = diag(p)';
end

%% Plot Linear Regression

%%% By Group %%%
figure; 
hold on;
for v = 1:2
    dataToPlot = squeeze(bLinGroup(v,subjects{v},:));
    plot(dataToPlot', 'o', 'color', color(v,:))
    errorbar(mean(dataToPlot),std(dataToPlot)./N(v),'color',color(v,:),'linewidth',2)
end
plot([0 5],[0 0],'k--') % plot reference line
set(gca,'xtick',1:length(linGroupRegressorNames),'xticklabels',linGroupRegressorNames)
ylabel('\beta  Weights')
%% Logistic Regression
% Try to predict item by item if they will get it correct 
regressorNames = {'Intercept','Span','First Item?' 'Last Item?', 'Equation Before Corr?'...
    'Span*EquationCorr'};
for v = 1:2 % loop through each group
    for su =subjects{v} % loop through each subject
        clear X Y rawData rawDataEqn
        %load this subject's raw data
        rawData = eval([vision{v},'Data(3).su',num2str(su),'.recall(:,9:10)']);
        rawDataEqn = eval([vision{v},'Data(3).su',num2str(su),'.eqtn(:,[10 14])']);
        rawDataEqn = cellfun(@(x) double(x),rawDataEqn); %convert to double
        numTrials = size(rawData,1);
        highestSpanReached = round(numTrials/3 + 1);
        % predictors: span, first item?, last item?, got equation after
        % correct?
        X = nan(sum(3.*span(1:highestSpanReached-1)), 5);
        Y = nan(sum(3.*span(1:highestSpanReached-1)),1);
        for trial = 1:size(rawData,1) % loop through trials
            % which letters did they correctly remember?
            compare = char(rawData(trial,:)); % the two we want to compare
            % delete any empty strings
            compare(:,compare(1,:) == ' ') = [];
            load = length(compare(1,:));
            eachLettCorrect = compare(1,:) == compare(2,:); % which letters were correctly remembered?
            % how did they preform on the equation before this
            % item
            eqnForThisTrial = rawDataEqn(rawDataEqn(:,2) == trial,1); % accuracy for the equtions of this trial
            leftTofill = find(isnan(Y));
            indToStartFill = leftTofill(1);
            if isempty(eqnForThisTrial) % did not collect some equation data, treat as if they missed it
                eqnForThisTrial = repmat([0],1,load);
            end
            for item = 1:load %loop through each item
                X(indToStartFill+(item-1),1:4) = [load (item == 1) (item==load)  eqnForThisTrial(item)];
            end
            Y(indToStartFill:indToStartFill+(load-1)) = eachLettCorrect;
        end
        % add interactions
        X(:,5) = [X(:,1).*X(:,4)]; 
        % z-score X
        X = zscore(X);
        % Must make Y have only positive values
        Y = Y+1;
        [bLog(v,su,:),~,~] = mnrfit(X,Y);
    end
end
%% Beta Weight Statistical Testing & Plotting

%%% statistical testing %%%

% Test against 0 for both groups
clear h p
for v = 1:2
    [h(v,:,:),p(v,:,:)]=ttest(bLog(v,subjects{v},:));
end

% Plot weights
figure
hold on;
for v = 1:2
    dataToPlot = squeeze(bLog(v,subjects{v},:))';
    plot(dataToPlot,'o','Color',color(v,:)) % plot beta weights for each subject
    errorbar(mean(dataToPlot'),std(dataToPlot')... % calculate thr within group mean & sd
    /sqrt(N(v)),'color',color(v,:),'linewidth',2)
end
plot([0 7], [0 0], 'k--') % reference line at y = 0
set(gca,'xtick',1:length(regressorNames),'xticklabels',regressorNames)
ylabel('\beta  Weights')
xlim([0 7])
%legend(vision)

% now test against each other
clear h p 
[h,p]=ttest2(bLog(1,subjects{1},:),bLog(2,subjects{2},:));

%% Group X Span X Eqtn Perf Anova 

% Group and Span are categorical variables. Equation performance is a continuous variable.
% Group is between subject and while span and eqtn perf are within subject (repeated measures) variables. 
% Trying to predict proportion correct positions averaged per span 
allSubsN = sum(N); 
allSubs = [subjects{1} subjects{2}];
spans= [2:9];
% Intialize Y and group vectors 
[Y,load,eqnPerf] = deal(nan(allSubsN*8,1)); % 8 data points for 42 subjects, 8*42 = 336
% deal function sets all the specified variables (your outputs) to your input argument
[group] = cell(allSubsN*8,1); 
% starting index for each subject for their 8 data points
startIndexPerSub = 1:8:allSubsN*8;
% Fill Y vector and group vectors
for sub = 1:allSubsN
    thisSubsInd = startIndexPerSub(sub): startIndexPerSub(sub)+7; % 8 indices where this subject's data lives
    isBlind = sub > N(1); % there are N(1) sighted subjects, subjects beyond this ind are blind
    group(thisSubsInd) = repmat({vision{isBlind+1}},1,8);
    for sp = 1:8 % loop through spans
        eqnPerf(thisSubsInd(sp)) = mathPerfperSpan(isBlind+1,allSubs(sub),spans(sp));
        load(thisSubsInd(sp)) = sp;
        Y(thisSubsInd(sp)) = eval(['intraSubAvg.',vision{isBlind+1},'(',...
            num2str(allSubs(sub)), ',1,',num2str(sp),');']);
    end
end

% Run the 3 way anova 
[p,tbl,stats] = anovan(Y,{group,eqnPerf,load},'model','full','continuous',[3],'varnames',{'group','load','eqnPerf'});