% WMSpan Project - Bedny Lab %%%%%%%%%%%%%%%
% 2018
% Data collected by Rita Loiotile
% Data analyzed by Rita Loiotile, Karen Arcos, and Nora Harhen
% Script by Nora Harhen (nharhen1@jhu.edu) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze across tasks. Run correaltions and do factor analysis. 
% Results are plotted and displayed to the command window. You can specifcy
% if you want the results displayed to be saved to a textfile by running the function
% setWriteToText(arg). arg = 1 write to text file. arg = 0 do not write to
% text file. 

% Load overall performance data from other tasks
clear all;
recog = load('recognitionOverallPerf.mat'); 
recall = load('recallOverallPerf.mat');
complex = load('complexSpanOverallPerf.mat'); 
clearvars -except recog recall complex
load writeToFile % specify whether or not you want outputs saved to text file using diary 

% Define subjects. Only including those who did all 6 tasks. 
subjects{1} = [1:10 12:16 18:23]; %S subjects
subjects{2} = [1:2 4:5 7:10 13:16 18:25]; %CB subjects
N = cellfun(@(x) length(x), subjects); % find number of subjects in each group
groups = {'S','CB'};
%% Set default settings for plotting 
% Because i am lazy and don't wanna set this every time
addpath('~/Desktop/Work/Research/WM/PlottingFunctions') % access the plotting functions I wrote
set(groot,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTitleFontSize',1)
set(groot,'DefaultLineMarkerSize',10)
set(groot,'DefaultLineLineWidth',2)
color = [244 140 66;66 155 244]./255; %sighted is second row (orange), cb is second row (blue)
%% Simple Recall and Complex Recall Correlations

%%% Forward and Backward Comparison %%%

% Calculate Correlations (Pearson) 
for group = 1:2
    for var = 1:2 % loop through proportion correct for positions and numbers
         toCorr = squeeze(eval(['recall.allSubOverallPerf.',groups{group},'(subjects{group},var,:)'])); 
        [rhoRecall(group,var),pRecall(group,var)] = corr(toCorr(:,1),toCorr(:,2)); 
    end
end

% Plot Results 
measures = {'% Positions Correct', '%  Letters Correct'};
simpleTasks = {'Forward Recall' ,' Backward Recall'};
figure
for plt = 1:2
    subplot(1,2,plt)
    hold on;
    % set up data
    sighted = squeeze(recall.allSubOverallPerf.S(subjects{1},plt,:)); 
    blind = squeeze(recall.allSubOverallPerf.CB(subjects{2},plt,:)); 
    % now plot 
    plotGroupCorrelation(sighted, blind, color, [0.4 1; 0.4 1], {simpleTasks{1} simpleTasks{2}}, measures{plt})
end

% Calculate difference between performance on forward and backward on our
% two measures then scatter plot the two groups overlaid on each other
diff = nan(2,25,2); % first dim: group, second: subjects, third: measure (positions or letters) 
diff(1,1:23,:)= recall.allSubOverallPerf.S(:,1:2,1) - recall.allSubOverallPerf.S(:,1:2,2); 
diff(2,1:25,:)= recall.allSubOverallPerf.CB(:,1:2,1) - recall.allSubOverallPerf.CB(:,1:2,2);

% Do a quick t-test on these 
[hDiff,pDiff,~,statDiff] = ttest2(diff(1,subjects{1},:), diff(2,subjects{2},:));

% Plot the results
figure
for measure = 1:2
    subplot(1,2,measure)
    hold on
    for group = 1:2
        plot(diff(group,subjects{group},measure),'o','Color',color(group,:))
    end
    ylabel('Forward Perf - Backward Perf')
    xlabel('Subjects')
    plot([0 30],[0 0], 'k--')
    title(measures{measure})
end

%%% Simple Recall and Complex Recall Comparison %%%

% Calculate Results
for group = 1:2 % sighted = 1 and blind = 2
    clear toCorr
    for task = 1:2 % forward = 1 and backward = 2
        for var = 1:2 % positions = 1 and letters = %2 
            toCorr(:,1)= eval(['recall.allSubOverallPerf.',groups{group},...
                '(subjects{group},var,task)']);
            toCorr(:,2)= eval(['complex.allSubOverallPerf.',groups{group},...
                '(subjects{group},var)']);
            [rhoMixRecall(group,task,var),pMixRecall(group,task,var)] = corr(toCorr(:,1), toCorr(:,2));
        end
    end
end

% Plot Results
figure
for task = 1:2
    for var = 1:2
        % set up data for plotting 
        sighted = [recall.allSubOverallPerf.S(subjects{1},var,task) ... %first col is simple recall, second col is complex
            complex.allSubOverallPerf.S(subjects{1},var)];
        blind =[recall.allSubOverallPerf.CB(subjects{2},var,task) ... % ditto'
            complex.allSubOverallPerf.CB(subjects{2},var)];
        % now plot 
        subplot(2,2,2*(var-1) + task)
        plotGroupCorrelation(sighted, blind, color, [0.2 1; 0.2 1], {simpleTasks{task} 'ComplexRecall'}, measures{var})
    end
end

%% Recognition Correlations
for group = 1:2
    toCorr = squeeze(eval(['recog.allSubOverallPerf.',groups{group},'(subjects{group},:)']));
    [rhoRecog(group),pRecog(group)] = corr(toCorr(:,1),toCorr(:,2));
end

figure
% set up data
sighted = squeeze(recog.allSubOverallPerf.S(subjects{1},:));
blind = squeeze(recog.allSubOverallPerf.CB(subjects{2},:));
% plot
plotGroupCorrelation(sighted, blind, color, [0.5 1; 0.6 1], {'Non-Verbal' 'Verbal'}, 'Recognition Tasks')

ylim([0.5 1]); xlim([0.6 1])
xlabel('Non-Verbal'); ylabel('Verbal')


%% Display these results to the command window and write to text file if wanted
if writeToFile
    diary('outputsAcrossTaskAnalysis')
end

%%% Backward and Forward Recall %%%
disp('These are the results of correlating (pearson correlation) forward recall performance with backward recall performance:')
disp('The first row is for sighted. The second is for blind. ')
disp('The first column is for % positions correct. The second is for % letters correct')
disp('Rho:')
disp(num2str(rhoRecall))
disp('P-values:')
disp(num2str(pRecall))

% Take the difference in performance for backward and forward for each
% participant. Do a t-test between groups for these differences. 
if sum(h) == 1
    disp('No signifigant difference in effect of manipulation for either positions correct or letters correct')
else
    disp('These was a signifigant difference of the effect of manipulation on these measures:')
    disp(measures(hDiff==1)) % display name of measure
    disp(num2str(pDiff(hDiff==1))) % display p-value 
    disp(num2str(statDiff.tstat(hDiff==1))) % display t-value 
end


%%% Simple and Complex Recall %%%
tasks = {'forward', 'backward'};
for t = 1:2
    disp(['These are the results of correlating (pearson correlation) simple ', tasks{t} ' recall' ...
        'performance with complex recall performance:'])
    disp('The first row is for sighted. The second is for blind. ')
    disp('The first column is for % positions correct. The second is for % letters correct')
    disp('Rho:')
    disp(num2str(rhoMixRecall(:,:,t)))
    disp('P-values:')
    disp(num2str(pMixRecall(:,:,t)))
end

%%% Recognition %%%
disp(['These are the results of correlating (pearson correlation) non-verbal recognition performance' ...
    'with verbal recognition performance:'])
disp('The first item is for sighted. The second is for blind.')
disp('Rho:')
disp(num2str(rhoRecog))
disp('P-values:')
disp(num2str(pRecog))

if writeToFile
    diary off 
end