% WMSpan Project - Bedny Lab %%%%%%%%%%%%%%%
% 2018
% Data collected by Rita Loiotile
% Data analyzed by Rita Loiotile, Karen Arcos, and Nora Harhen
% Script by Nora Harhen (nharhen1@jhu.edu)  and Karen Arcos
% (karcos1@uci.edu) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze Nonverbal discern pairs task data for ensuring that sighted and 
% CB participants are matched on auditory discernment as a control for the
% recognition tasks. 

clear all;
load('allData_21-Aug-2018 15:08:37.mat')
%% Load the data and combine

%%% Initialize Important Variables %%%
vision = {'S', 'CB'};  % our two groups
subjects{1} = [1:10 12:23]; %S subjects, eliminated o w/o critical tasks
subjects{2} = [1:2 4:5 7:10 13:16 18:22 24:25]; %CB subjects, eliminated ones w/o critical tasks
N = cellfun(@(x) length(x), subjects); % find number of subjects in each group

toKeep = [10,13]; % columns we care about: isMatchTarget and isAnswerCorrect

overallPerf = nan(25,2); % intialize 

%%% Commence the looping and the Combining! %%%

for v = 1:length(subjects) % loop through our two groups 
    for su = subjects{v} % loop through the subjects within the current group
            tempMat(:,1) = eval(['cell2mat(',vision{v},'Data(',num2str(4),').su',num2str(su),'(:,toKeep(1)))']);
            tempMat(:,2) = eval(['cell2mat(',vision{v},'Data(',num2str(4),').su',num2str(su),'(:,toKeep(2)))']);

            % get overall performance measure
            overallPerf(su,v) = mean(tempMat(:,2));
             
            % get sensitivity measures
            criteraConfusionMat = [0 0; 0 1; 1 1; 1 0]; % false positives,correct misses, hits, false negatives
            [tf, index] = ismember(tempMat, criteraConfusionMat, 'rows');  % find them
            conMat = arrayfun(@(x) sum(ismember(index,x)),1:4);  % count them up for confusionMatrix
            conMat = conMat./sum(conMat); % get a frequency 
            conMat = reshape(conMat,2,2); % reshape from vector to matrix 
            eval(['sensitivites.', vision{v}, '(:,:,su)= conMat;']) % save data 
    end
end

%% Calculate Group stats

% calculate stats, then plot and display data 
statNames = {'mean','std','min', 'max','ci'}; 

% Set up variables for calculating confidence interval
T_mult = abs(tinv([0.025 0.975], N-1));%find the t-stat multiplier, dependent on group size

for s = 1:length(statNames)-1 % loop through stats to compute, leave ci out
    eval(['stat(s,:)= [', ... 
        statNames{s},'(overallPerf(subjects{1},1)) ', statNames{s},'(overallPerf(subjects{2},2))];'])
end % for s

% do ci
stat(5,:) = T_mult(v)*stat(2,:)./sqrt(N);

%% Set default settings for plotting 
% Because i am lazy and don't wanna set this every time
set(groot,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTitleFontSize',1)
color = [244 140 66;66 119 244]./255; %s is first row (orange) , cb is second row (blue) 
%% Plot the performance %%%
figure
hold on;
for b = 1:2
    bar(b,stat(1,b),'facecolor', color(b,:))
    errorbar(b,stat(1,b),stat(2,b),'.k')
    set(gca,'xtick',1:2,'xticklabels',{'S','CB'})
end
xlabel('Group'); ylabel('% Correct')

%%% Statistical Testing %%%

[h,p,~,stats] = ttest2(overallPerf(subjects{1},1), overallPerf(subjects{2},2)); 
