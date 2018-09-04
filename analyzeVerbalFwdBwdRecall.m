%%% Analyze Verbal Forward Recall and Verbal Backward Recall %%%%

% Script by: Nora Harhen (nharhen1@jhu.edu) 
% Bedny Lab
% August 2018
clear all;
load('allData_28-Aug-2018 11:27:13.mat')
%% Collapse all the data into two matrices, one for each task

% Our two groups 
vision = {'S', 'CB'}; 
span = [2:9];

subjects{1} = [1:10 12:16 18:23]; %S subjects, eliminated ones w/o data and w/o critical tasks
subjects{2} = [1:2 4:5 7:10 13:16 18:25]; %CB subjects, eliminated ones w/o data and w/o critical tasks
N = cellfun(@(x) length(x), subjects); % find number of subjects in each group 

% Columns we care about for this analysis
toKeep = [13,15:16]; % proportion correct positions, prop corr letters, rxn time
toKeepNonNumer = [9:10];
headingsToKeep = cell2mat(headings{1}(toKeep)); %headings are same for both tasks

% initialize matrix 

% Create matrix of data
for v = 1:length(subjects) % loop through our two groups
    for su = subjects{v} %loop through subjects in groups
        for task = 1:2 % loop through tasks
            % create temporary matrix 
            tempMat = eval(['cell2mat(',vision{v},'Data(',num2str(task),').su',num2str(su),'(:,toKeep))']);
            if size(tempMat,1) < 16 % need to zero pad if reached less than span 9
                numRowsAdd = 16 - size(tempMat,1);
                rowsToAdd= [0 0 -1];
                tempMat = [tempMat;repmat(rowsToAdd,numRowsAdd,1)];
            end
            % average across span, each span has two trials associated with
            % it, so we can average every other row
            n=2;
            averages = arrayfun(@(i) mean(tempMat(i:i+n-1,:)),... % averages is a cell array 
                1:n:length(tempMat)-n+1,'UniformOutput', false);
            averages = cell2mat(averages); % convert from cell to matrix 
            % reshape data from vector to 3x8 matrix 
            % save to giant data matrix 
            eval(['intraSubAvg.', vision{v},'(',num2str(su),',:,:,', num2str(task),') = reshape(averages,[3,8]);']);
            
        end
    end
end
%% Set default settings for plotting 
% Because i am lazy and don't wanna set this every time
set(groot,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTitleFontSize',1)
set(groot,'DefaultLineMarkerSize',10)
set(groot,'DefaultLineLineWidth',2)
color = [244 140 66;66 155 244]./255; %sighted is second row (orange), cb is second row (blue)
%% Average across subjects

%%% Performance Across Spans  %%%
[statPerSpan(1)] = descStats(intraSubAvg.S(subjects{1},1:2,:,:));
[statPerSpan(2)] = descStats(intraSubAvg.CB(subjects{2},1:2,:,:));

figure;
titles = {'Forward Recall', 'Backward Recall'};
yaxisLabels = {'Proportion Positions Correct (%)', 'Proportion Letters Correct (%)'};
for i = 1:2 % Loop through positions versus letters correct
    for j = 1:2 % loop through forward and backward recall
        subplot(2,2,i+j*(j-1))
        hold on;
        for v = 1:2
            points = squeeze(statPerSpan(v).mean(:,i,:,j)); % markers to plot
            errBars = squeeze(statPerSpan(v).ci(:,i,:,j)); % error bars to plot
            errorbar(points,errBars,'o', 'MarkerFaceColor',color(v,:), 'Color', color(v,:))
        end
        ylabel(yaxisLabels{i}); xlabel('Span Length')
        set(gca,'xtick',1:8,'xticklabels',2:9)
        title(titles{j})
        ylim([0 1])
    end
    legend('S', 'CB')
end

%%%  Statistical testing per each span and each task %%%
[h,p,stats.tstat] = ttest2(intraSubAvg.CB(subjects{2},:,:,:), intraSubAvg.S(subjects{1},:,:,:));
for task = 1:2 % Note: for backward, no one got to span 9 so they are all 0, which is why it is nan in the last column
    disp(['These are the p-values for ', titles{task},'.'])
    disp('The first row is for the proportion of positions correct. The second is for the proportion of letters correct.')
    squeeze(p(:,1:2,:,task))%forward
end

%%% Plot overall performance %%%

% Average across spans per subject
overallPerfPerSub.S = squeeze(mean(intraSubAvg.S(subjects{1},:,:,:),3)); % new size: 21 x 3 x 2
overallPerfPerSub.CB = squeeze(mean(intraSubAvg.CB(subjects{2},:,:,:),3)); % new size: 21 x 3 x 2

% Average across subjects within group
[statOverall(1)] = descStats(overallPerfPerSub.S(:,1:2,:));
[statOverall(2)] = descStats(overallPerfPerSub.CB(:,1:2,:));

% Overlay the plots for CB and S showing performance for forward and
% backward recall for proportion positions correct and proportion letters
% correct
lineStyles = {'o-', 'o--'}; % solid line for positions and dashed for letters 

figure;
hold on;
for v = 1:2 %lop through group
    for measure = 1:2 % loop through measures of performance (positions and letters)
        % select data to plot
        marker = squeeze(statOverall(v).mean(:,measure,:)); 
        errBar = squeeze(statOverall(v).ci(:,measure,:));
        % plot data
        e(v) = errorbar(marker,errBar,lineStyles{measure}); 
       [e(v).Color, e(v).MarkerFaceColor] = deal(color(v,:));

    end
end
ylabel('Proportion correct');
xlim([0 3]); ylim([0.45 1])
set(gca,'xtick',1:2,'xticklabels',{'Forward', 'Backward'})
set(gca, 'fontsize',16)
legend('S Positions', 'S Letters','CB Positions', 'CB Letters')
title('Overall Performance')


%%% statistical testing %%%
[h,p,stats.tstat] = ttest2(overallPerfPerSub.CB, overallPerfPerSub.S);
for task = 1:2
    disp(['These are the p-values from a t-test between groups for ', titles{task},' overall performance.'])
    disp('The first column is for positions and the second is for letters')
    p(1,1:2,task)
    disp(['These are the t-stats from a t-test between groups for ', titles{task},' overall performance.'])
    disp('The first column is for positions and the second is for letters')
    stats.tstat(1,1:2,task)
end

%%% Plot Effect Sizes %%%

% Effect of Manipulation (Forward vs. Backward) 
effectManip.CB = overallPerfPerSub.CB(:,1:2,1) - overallPerfPerSub.CB(:,1:2,2);
effectManip.S = overallPerfPerSub.S(:,1:2,1) - overallPerfPerSub.S(:,1:2,2);

[effectManipOverall(1)] = descStats(effectManip.S);
[effectManipOverall(2)] = descStats(effectManip.CB);

titles={'Positions', 'Letters'}; 

figure;
for plt = 1:2
    subplot(1,2,plt)
    hold on;
    for v =1:2
        dataToPlot = effectManipOverall(v).mean(plt); % the bar
        intervals = effectManipOverall(v).ci(plt); % confidence intervals
        bar(v,dataToPlot,'FaceColor',color(v,:))
        errorbar(v,dataToPlot,intervals, '.k')
    end
    xlim([0 3]); ylim([0 0.3])
    set(gca,'xtick',1:2,'xticklabels',{'CB', 'S'})
    title(titles{plt});ylabel('Effect of Manip: Perf Fwd - Perf Bwd')
end


% Effect of Manipulation over all Spans 

% Find the difference in performance
effectManipAllSpan.S = intraSubAvg.S(subjects{1},1:2,:,1) - intraSubAvg.S(subjects{1},1:2,:,2); 
effectManipAllSpan.CB = intraSubAvg.CB(subjects{2},1:2,:,1) - intraSubAvg.CB(subjects{2},1:2,:,2); 

% Average the differences across participants 
[effectManipPerSpan(1)] = descStats(effectManipAllSpan.S );
[effectManipPerSpan(2)] = descStats(effectManipAllSpan.CB);

%%% Plot the data %%%
figure;
for plt = 1:2
    subplot(1,2,plt)
    points = [effectManipPerSpan(1).mean(1,plt,:) ;effectManipPerSpan(2).mean(1,plt,:)];
    errBars = [effectManipPerSpan(1).ci(1,plt,:) ;effectManipPerSpan(2).ci(1,plt,:)];
    
    hBar=bar(squeeze(points)',1);
    xBar=cell2mat(get(hBar,'XData')).' + [hBar.XOffset]; % compute bar centers
    hold on;
    hEB=errorbar(xBar,squeeze(points)',squeeze(errBars)','.k');
    for i = 1:2 % set colors by hand
        set(hBar(i),'FaceColor',color(i,:))
    end
    title(titles{plt});ylabel('Effect of Manip: Perf Fwd - Perf Bwd')
    set(gca,'xtick',1:8,'xticklabels',2:9)
    legend('S','CB')
    xlabel('Span')
end

%%% Statistical testing %%%
for v =1:2
    strInput = ['effectManipAllSpan.', vision{v}];
    [h,p,stats.tstat]= ttest(eval(strInput)); %within group ttest, compare to 0 
    disp(['T-testing the effect of manipulation in the ', vision{v}, 'group against 0:'])
    disp('These are the p-values per spans 2-9. First row is positions and second is letters')
    squeeze(p)
end

[h,p,stats.tstat]= ttest2(effectManipAllSpan.CB,effectManipAllSpan.S); % between group ttest 
disp(['T-testing the effect of manipulation between the two groups:'])
disp('These are the p-values per spans 2-9. First row is positions and second is letters')
squeeze(p)

%% Group X Task X Span Anova 

% Group, Task, and Span are categorical variables. Group is between subject and
% Task and Span are within subject (repeated measures) variables. 
allSubsN = sum(N); 
allSubs = [subjects{1} subjects{2}];
taskNames = {'Forward', 'Backward'}; 
% Intialize Y and group vectors 
[Y,load,subNo] = deal(nan(allSubsN*16,1)); % 16 data points for 41 subjects
% deal function sets all the specified variables (your outputs) to your input argument
[group,task] = deal(cell(allSubsN*16,1)); 
% starting index for each subject for their 16 data points
startIndexPerSub = 1:16:allSubsN*16;
% Fill Y vector and group vectors
for sub = 1:allSubsN
    thisSubsInd = startIndexPerSub(sub): startIndexPerSub(sub)+15; % 8 indices where this subject's data lives
    isBlind = sub > N(1); % there are N(1) sighted subjects, subjects beyond this ind are blind
    group(thisSubsInd) = repmat({vision{isBlind+1}},1,16); 
    %subNo(thisSubsInd) = repmat(sub,1,8); 
    for t = 1:2 % loop through task
        currTaskInd = thisSubsInd(8*t-7:8*t); 
        task(currTaskInd) = repmat({taskNames{t}},1,8);
        for sp = 1:8 % loop through spans
            load(currTaskInd(sp)) = sp; 
            Y(currTaskInd(sp)) = eval(['intraSubAvg.',vision{isBlind+1},'(',...
                num2str(allSubs(sub)), ',1,',num2str(sp),',',num2str(t),');']);
        end
    end
end

% Run the 3 way anova 
[p,tbl,stats] = anovan(Y,{group,task,load},'model','full','varnames',{'group','task','span'});