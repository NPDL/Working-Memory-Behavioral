% WMSpan Project - Bedny Lab %%%%%%%%%%%%%%%
% 2018
% Data collected by Rita Loiotile
% Data analyzed by Rita Loiotile, Karen Arcos, and Nora Harhen
% Script by Nora Harhen (nharhen1@jhu.edu)  and Karen Arcos
% (karcos1@uci.edu) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze Recognition task data for verbal and nonverbal as well as data
% from nonverbal discern pairs task as a control 

clear all;
addpath('~/Desktop/Work/Research/WM') % access files outside this subfolder
load('allData_28-Aug-2018 11:27:13.mat')
load('writeToFile') % specifies if you want what is output to command window to be saved to text file
if writeToFile
    diary('outputsRecognitionAnalysis')
end
%% Clean/reshape data and merge together

% A short history of this data set and what needs to be corrected

% Data was collected differently from 2014 & 2016 to 2018. During 2014 &
% 2016, participants were stopped following a span in which their
% performance was below 50%. The spans went from 3 to 15 for nonverbal and
% 5 to 15 for verbal. In 2018, the spans were trimed to 3 to 6 for
% nonverbal and 5 to 8 for verbal.

%%% Initialize Important Variables %%%

subjects{1} = [1:10 12:23]; %S subjects, eliminated o w/o critical tasks
subjects{2} = [1:2 4:5 7:10 13:16 18:25]; %CB subjects, eliminated ones w/o critical tasks
N = cellfun(@(x) length(x), subjects); % find number of subjects in each group

toKeepNumer = [10,13]; % columns of numeric data we want to keep from the raw data files
%  13 = correct/not correct answer
toKeepNonNumer = [11]; % same as above but with non-numeric data, lure type per trial
vision = {'S', 'CB'};  % our two groups
lureTypes  = {'IdentityChange' ,'N/A','SlideOneOver','SwapTwo', 'NotReached'};
taskNames = {'Nonverbal', 'Verbal'};
% Track number of partipants who's span extended beyond or below our bounds
% in a 2 x 2 X 6 matrix: first dim: num below or above (below is first row, above is second row),
% second dim: group from previous years, third dim: task
spanOutBounds = zeros(2,2,6);

%%% Loop through groups, subjects, and tasks to compile data %%%

% should be 32 trials total
for v = 1:length(subjects) % loop through our two groups
    for su = subjects{v} % loop through the subjects within the current group
        for task = 5:6 % loop through our two tasks: nonverbal & verbal recognition
            clear tempMat cellMatNonNumer
            % first make a cell array of the columns we want
            tempMat(:,1) = eval(['cell2mat(',vision{v},'Data(',num2str(task),').su',num2str(su),'(:,toKeepNumer(1)))']);
            tempMat(:,2) = eval(['cell2mat(',vision{v},'Data(',num2str(task),').su',num2str(su),'(:,toKeepNumer(2)))']);
            cellMatNonNumer = eval([vision{v},'Data(',num2str(task),').su',num2str(su),'(:,toKeepNonNumer)']);
            % then convert to a matrix of type double
            
            if size(tempMat,1) < 32 % need to zero pad if reached less than span 6/8
                disp([vision{v},num2str(su), ' needed to zero pad'])
                % get sensitivity measures before zero padding
                criteraConfusionMat = [1 1;1 0;0 0;0 1]; % hits, false negatives, false alarms, correct rejection
                [tf, index] = ismember(tempMat, criteraConfusionMat, 'rows');  % find them
                conMat = arrayfun(@(x) sum(ismember(index,x)),1:4);  % count them up for confusionMatrix
                nTrueMatches = sum(conMat(1:2)); nTrueNonMatches = sum(conMat(3:4));
                conMat = conMat./[nTrueMatches nTrueMatches nTrueNonMatches nTrueNonMatches]; % get a frequency
                conMat = reshape(conMat,2,2); % reshape from vector to matrix
                eval(['sensitivites.', vision{v}, '(:,:,su,task)= conMat;']) % save data
                
                % Numeric Data
                numRowsAdd = 32 - length(tempMat);
                rowsToAdd= [ones(numRowsAdd,1).*[-1 0.5]];
                tempMat = [tempMat; rowsToAdd];
                
                % Non-numeric Data
                cellsToAdd = cell(1,numRowsAdd);
                cellsToAdd(:,:) = {'NotReached'};
                cellMatNonNumer = [cellMatNonNumer ;cellsToAdd'];
                
                % update counts
                spanOutBounds(1,v,task) = spanOutBounds(1,v,task) + 1;
            elseif size(tempMat,1) > 32 % need to trim the data if went beyond span 6/8
                disp([vision{v},num2str(su), ' needed to trim'])
                % get sensitivity measures before trimming
                criteraConfusionMat = [1 1;1 0;0 0;0 1]; % hits, false negatives, false alarms, correct rejection
                [tf, index] = ismember(tempMat, criteraConfusionMat, 'rows');  % find them
                conMat = arrayfun(@(x) sum(ismember(index,x)),1:4);  % count them up for confusionMatrix
                nTrueMatches = sum(conMat(1:2)); nTrueNonMatches = sum(conMat(3:4));
                conMat = conMat./[nTrueMatches nTrueMatches nTrueNonMatches nTrueNonMatches]; % get a frequency, divide false alarms and hits by (false alarms + hits) and
                conMat = reshape(conMat,2,2); % reshape from vector to matrix
                eval(['sensitivites.', vision{v}, '(:,:,su,task)= conMat;']) % save data
                
                tempMat = tempMat(1:32,:);
                cellMatNonNumer = cellMatNonNumer(1:32);
                % update counts
                spanOutBounds(2,v,task) = spanOutBounds(2,v,task) + 1;
            elseif size(tempMat,1) == 32
                disp([vision{v},num2str(su), ' perfect size'])
                % get sensitivity measures before zero padding
                criteraConfusionMat = [1 1;1 0;0 0;0 1]; % hits, false negatives, false alarms, correct rejection
                [tf, index] = ismember(tempMat, criteraConfusionMat, 'rows');  % find them
                conMat = arrayfun(@(x) sum(ismember(index,x)),1:4);  % count them up for confusionMatrix
                nTrueMatches = sum(conMat(1:2)); nTrueNonMatches = sum(conMat(3:4));
                conMat = conMat./[nTrueMatches nTrueMatches nTrueNonMatches nTrueNonMatches]; % get a frequency, divide false alarms and hits by (false alarms + hits) and
                conMat = reshape(conMat,2,2); % reshape from vector to matrix
                eval(['sensitivites.', vision{v}, '(:,:,su,task)= conMat;']) % save data
            end
            % Now average across each span
            n=8;
            averages = arrayfun(@(i) mean(tempMat(i:i+n-1,2)),... % averages is a cell array
                1:n:length(tempMat)-n+1,'UniformOutput', false);
            % convert to matrix from cell array
            averages = cell2mat(averages);
            % if performance is below 0.5 for a certain span then set to
            % 0.5
            averages(averages < 0.5) = 0.5;
            %save to giant data matrix that will be a part of a structure
            eval(['intraSubAvg.', vision{v},'(',num2str(su),',:,', num2str(task),') = averages;']);
            
            % select the error trials for lure counts
            luresWorked = cellMatNonNumer(tempMat(:,2)<1); % Participants made mistake
            % Get lure counts and save to matrix which is part of a
            % structure that has fields for each group
            getCounts = cellfun(@(x) sum(ismember(luresWorked,x)),lureTypes,'un',0);
            getCounts=cell2mat(getCounts); % convert to matrix
            eval(['lureCounts.', vision{v},'(',num2str(su),',:,', num2str(task),') = getCounts;']);
        end
    end
end
% delete unnecessary empty space
intraSubAvg.S = intraSubAvg.S(:,:,5:6);
intraSubAvg.CB = intraSubAvg.CB(:,:,5:6);
sensitivites.S=sensitivites.S(:,:,:,5:6);
sensitivites.CB=sensitivites.CB(:,:,:,5:6);
%% %% Set default settings for plotting 
% Because i am lazy and don't wanna set this every time
set(groot,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesTitleFontSize',1)
set(groot,'DefaultLineMarkerSize',10)
set(groot,'DefaultLineLineWidth',2)
color = [244 140 66;66 119 244]./255; %s is first row (orange) , cb is second row (blue) 
addpath('~/Desktop/Work/Research/WM/PlottingFunctions') % access the plotting functions I wrote
%% Calculate averages across subjects (within groups), standard deviation, and ci 

%%% Performance Across Spans %%%

% Calc Stats %
[statPerSpan(1)] = descStats(intraSubAvg.S(subjects{1},:,:));
[statPerSpan(2)] = descStats(intraSubAvg.CB(subjects{2},:,:));

% Display to the command window 
disp('Stats for performance separated by span')
disp('Span increases from left to right - Non-Verbal: 3-6, Verbal Spans: 5-8')
fieldNames = fieldnames(statPerSpan);
for task=1:2
    disp(taskNames{task})
    for s = 1:length(fieldNames)
        disp(fieldNames{s})
        for group = 1:2
            disp(vision{group})
            dataToDisplay = getfield(statPerSpan,{group},fieldNames{s},{1,1:4,task});
            disp(num2str(dataToDisplay))
        end
    end
end

% Plot Data %
spans(1,:) = 3:6; spans(2,:) = 5:8; 
titles = {'Non Verbal Fwd Recognition', 'Verbal Fwd Recognition'};

figure;
hold on
for plt = 1:2 %loop through the two tasks
    subplot(1,2,plt);
    hold on;
    for v = 1:2 % loop through the two groups 
         e(v) = errorbar(statPerSpan(v).mean(:,:,plt),statPerSpan(v).ci(:,:,plt), 'o'); % plot the marker and errorbar
         e(v).Color = color(v,:); % set color of errorbars 
         e(v).MarkerFaceColor = color(v,:);  % set color of markers
         e(v).MarkerSize = 10;  % set size of markers 
    end % for v
    set(gca,'xtick',1:4,'xticklabels',spans(plt,:)) % list the spans 
    ylabel('% Correct'); xlabel('Span'); title(titles{plt}); 
    legend('S','CB')
end % for plt

% Statistical Testing %

% T-Test for nonverbal and verbal
[h,p,~,stats] = ttest2(intraSubAvg.S(subjects{1},:,:), intraSubAvg.CB(subjects{2},:,:));

disp('Do a T-Test to compare verbal and non-verbal across all spans:')
for t = 1:2 % loop through tasks
    sigSpan = h(:,:,t)==1; 
    disp(titles{t})
    disp(['Statistically Signifigant Spans: ', num2str(spans(t,find(sigSpan)))])  
    disp(['p: ', num2str(p(:,sigSpan,t))])
    disp(['t-stat: ', num2str(stats.tstat(:,sigSpan,t))])
end

%%% Overall Performance %%%

% Calc Stats %
[statOverall(1)] = descStats(squeeze(mean(intraSubAvg.S(subjects{1},:,:),2)));
[statOverall(2)] = descStats(squeeze(mean(intraSubAvg.CB(subjects{2},:,:),2)));

% Display the Results 
disp('These are the descriptive stats for overall performance analysis')
fieldNames = fieldnames(statOverall);
for task = 1:2 % nonverbal = 1 verbal = 2
    disp(taskNames{task})
    for s = 1:length(fieldNames) % loop through fields which are different tasks
        disp(fieldNames{s})
        for group = 1:2 % sighted = 1 blind = 2
            dataToDisplay(group) = getfield(statOverall,{group},fieldNames{s},{task}); % select the struct at index group of the struct array statOverll
            % then select the fieldname specified by fieldNames{s}, select
            % the index task from the field called fieldNames{s}
        end
        disp(vision)
        disp(num2str(dataToDisplay))
    end
end

% Plot the Data %
figure;
hold on;
for v = 1:2 %loop through group
    % select data to plot
    marker = squeeze(statOverall(v).mean);
    errBar = squeeze(statOverall(v).ci);
    % plot data
    e(v)= errorbar(marker,errBar,'o-', 'MarkerFaceColor',color(v,:));
    [e(v).Color, e(v).MarkerFaceColor] = deal(color(v,:));
end
ylabel('% correct');
xlim([0.5 2.5]); ylim([0.5 1])
set(gca,'xtick',1:2,'xticklabels',{'NonVerbal', 'Verbal'})
set(gca, 'fontsize',16)
legend('S', 'CB')
title('Overall Performance')

% Statistical Testing % 
[h,p,~,stats] = ttest2(squeeze(mean(intraSubAvg.S(subjects{1},:,:),2)),squeeze(mean(intraSubAvg.CB(subjects{2},:,:),2)));
disp('Did t-test between groups on overall performance. Found signifigant differences on :')
disp(titles(find(h==1)))
disp(['P-value: ', num2str(p(h==1)), ' T-value: ', num2str(stats.tstat(h==1))])
%% Check differences in frequency of type of mistakes made
% see if the two groups tend to make the same mistakes or different
[statMistake(1)] = descStats(lureCounts.S(subjects{1},:,5:6));
[statMistake(2)] = descStats(lureCounts.CB(subjects{2},:,5:6));

% Display the Results 
disp(' These are the descriptive stats for the analysis of the mistakes')
disp('The mistakes are counts no proportions')
disp('The first row are the sighted. The second is the blind')
fieldNames = fieldnames(statMistake(1));
clear dataToDisplay
for task = 1:2 % nonverbal = 1 verbal = 2
    disp(taskNames{task})
    for s = 1:length(fieldNames) % loop through fields which are different tasks
        disp(fieldNames{s})
        for group = 1:2 % sighted = 1 blind = 2
            dataToDisplay(group,:) = getfield(statMistake,{group},fieldNames{s},{1,1:5,task}); % select the struct at index group of the struct array statOverll
            % then select the fieldname specified by fieldNames{s}, select
            % the index task from the field called fieldNames{s}
        end
        disp(lureTypes)
        disp(num2str(dataToDisplay))
    end
end

figure;
for plt = 1:2 % loop through tasks to plot
    subplot(1,2,plt)
    points = [statMistake(1).mean ;statMistake(2).mean];
    errBar = [statMistake(1).ci ;statMistake(2).ci];
    
    hBar=bar(points(:,:,plt)',1);
    xBar=cell2mat(get(hBar,'XData')).' + [hBar.XOffset]; % compute bar centers
    hold on
    hEB=errorbar(xBar,points(:,:,plt)',errBar(:,:,plt)','.k');
    set(gca,'xtick',1:5,'xticklabels',lureTypes)
    title(titles{plt}); ylabel('Number of Mistakes Made')
    for i = 1:2 % set colors by hand
        set(hBar(i),'FaceColor',color(i,:))
    end
legend('S','CB')
end

%%% Statistical Testing %%%
[h,p,~,stats] = ttest2(lureCounts.S(subjects{1},:,5:6),lureCounts.CB(subjects{2},:,5:6)); 
for task = 1:2
    disp(['Doing a t-test on ' taskNames{task},': LureCounts'])
    disp('The statistics for each type of lure are presented in this order:')
    disp(lureTypes)
    disp(['P-values: ', num2str(p(:,:,task))])
    disp(['T-values: ', num2str(stats.tstat(:,:,task))])
end

%%% Oversensitive vs. Undersensitive %%% 

% calculate rates for confusion matrix 
avgSense(:,:,:,1) = squeeze(mean(sensitivites.S(:,:,subjects{1},:),3));
avgSense(:,:,:,2) = squeeze(mean(sensitivites.CB(:,:,subjects{2},:),3));

% display the results
disp('These are the matrices for hits, false negatives, false alarams, and correct rejection')
disp('The first row are hits and then false negatives')
disp('The second row are false alarms and then correct rejections')
for task = 1:2
    disp(taskNames{task})
    for group = 1:2
        disp(vision{group})
        disp(num2str(avgSense(:,:,task,group)))
    end
end

clear h p
[h,p,~,stat] = ttest2(sensitivites.S(1,1,subjects{1},2),sensitivites.CB(1,1,subjects{2},2));
if h == 1
    disp('There is a signifigant difference in hit rates')
else 
    disp('There is no signifigant difference in hit rates')
end
disp(['Pvalue: ', num2str(p)]) 
disp(['Tvalue: ', num2str(stat.tstat)])

xlabels={'True Match' ,'True NonMatch'}; ylabels={'Chose Match' ,'Chose NonMatch'};
figure
for t = 1:2 % loop through task
    for g = 1:2 % loop through group
        subplot(2,2,t*(t-1) + g)
        heatmap(xlabels,ylabels,avgSense(:,:,t,g),'colormap',parula,'ColorLimits',[0 .5]);
        title([titles{t}, ': ', vision{g}])
    end
end

%% Calculate d' for both groups 

%initialize d prime mat to fill
dPrime = nan(2,25,2); % Group X Subjects in Group X Task

% Calculate dPrime
for v =1:2 % loop through group
    for t = 1:2 % loop through task
        hRate = eval(['squeeze(sensitivites.', vision{v},'(1,1,subjects{v},t));']); % hits
        faRate = eval(['squeeze(sensitivites.', vision{v},'(1,2,subjects{v},t));']);% false alarms
        % need to find  p = 0 or 1 and replace with 1/N or N-1/N, N (number of trials) = 32
        hRate(hRate == 0) = 1/32;
        hRate(hRate == 1) = 31/32;
        faRate(faRate == 0) = 1/32;
        faRate(faRate == 1) = 1/32;
        % z-transform and calculate d'
        dPrime(v,subjects{v},t) = norminv(hRate) - norminv(faRate);
    end
end

% Difference in DPrime between tasks
effectStimDPrime = nan(25,2); 
for v =1:2
    effectStimDPrime(subjects{v},v) = dPrime(v,subjects{v},1) - dPrime(v,subjects{v},2);
end

%%% Statistical Testing %%%
[h,p,~,stats] = ttest2(dPrime(1,subjects{1},:),dPrime(2,subjects{2},:));
disp('There is a signifigant difference between groups on dPrime for these tasks:')
disp(taskNames(h==1))
disp(['P-value: ', num2str(p(h==1)), ' T-value: ', num2str(stats.tstat(h==1))])
[h,p,~,stats] = ttest2(effectStimDPrime(subjects{1},1),effectStimDPrime(subjects{2},1));
if h == 1
    disp('There is a signifigant difference between groups on the effect of task on dPrime ')
    disp(['P-value: ', num2str(p(h==1)), ' T-value: ', num2str(stats.tstat(h==1))])
else
    disp('There is no signifigant difference between groups on the effect of task on dPrime ')
    disp(['P-value: ', num2str(p(h==0)), ' T-value: ', num2str(stats.tstat(h==0))])

end

%% Spans above and below our testing 
% Frequency particpants who were above or below the highest tested span in
% previous years 
NOld = [20 11]; % counted from spreadsheet 
spanOutBounds = spanOutBounds(:,:,5:6)./NOld;
%% Plot the Out of Bound Span Frequencies per group 
% Plot Something that looks like a confusion matrix 
ylabels = {'Below Highest Tested Span', 'Above Highest Tested Span'};
figure;
for plt = 1:2
    subplot(1,2,plt)
    h = heatmap(vision, ylabels, spanOutBounds(:,:,plt)); 
    h.Colormap = parula;
    h.ColorLimits = [0 1]; 
    title(titles{plt})
end

%% Do Linear Regression 
% Regressors: isBlind, span, typeStimuli, isBlind*span,
% isBlind*typeStimuli,isBlind*span*typeStimuli
% Predicting % Correct Choices 

linRegressorNames = {'Intercept', 'isBlind', 'Span', 'Type Stimuli', 'Match/Non-Match Trial', 'isBlind*Span',...
    'isBlind*typeStimuli','isBlind*Match/Non-Match','isBlind*Span*typeStimuli'};

% Combine data from each task for each group for each subject
allSubsN = sum(N); 
allSubs = [subjects{1} subjects{2}];

% for su = 1:allSubsN % loop through each subject
%     isBlind = (su > 22);
%     for task = 1:2 % loop through non-verbal and verbal
%         for load = 1:4 % four different loads tested in each tasks
%             
%         end
%     end
%     
% end

%% Group X Task X Span Anova 

% Group, Task, and Span are categorical variables. Group is between subject and
% Task and Span are within subject (repeated measures) variables. 
allSubsN = sum(N); 
allSubs = [subjects{1} subjects{2}];
% Intialize Y and group vectors 
[Y,span,subNo] = deal(nan(allSubsN*8,1)); % 8 data points for 42 subjects, 8*42 = 336
% deal function sets all the specified variables (your outputs) to your input argument
[group,task] = deal(cell(allSubsN*8,1)); 
% starting index for each subject for their 8 data points
startIndexPerSub = 1:8:allSubsN*8;
% Fill Y vector and group vectors
for sub = 1:allSubsN
    thisSubsInd = startIndexPerSub(sub): startIndexPerSub(sub)+7; % 8 indices where this subject's data lives
    isBlind = sub > N(1); % there are N(1) sighted subjects, subjects beyond this ind are blind
    group(thisSubsInd) = repmat({vision{isBlind+1}},1,8); 
    subNo(thisSubsInd) = repmat(sub,1,8); 
    for t = 1:2 % loop through task
        currTaskInd = thisSubsInd(4*t-3:4*t); 
        task(currTaskInd) = repmat({taskNames{t}},1,4);
        for load = 1:4 % loop through spans
            span(currTaskInd(load)) = load; 
            Y(currTaskInd(load)) = eval(['intraSubAvg.',vision{isBlind+1},'(',...
                num2str(allSubs(sub)), ',',num2str(load), ',',num2str(t),');']);
        end
    end
end

% Run the 3 way anova 
[p,tbl,stats] = anovan(Y,{group,task,span},'model','full','varnames',{'group','task','span'});
disp('Results of Anova: Group(2) X Task(2) X Span(3)')
disp(tbl)
%% Do correlation between discern pair performance and overall performance on each recognition task 
clear load p % clear variable so can use fucntion
% load discern pair data  for correlation 
load('discernPairsOverallPerf.mat')

% do separate correlations and plotting per group
for task = 1:2
    [rho(1,task),p(1,task)] = corr(squeeze(mean(intraSubAvg.S(subjects{1},:,task),2)), overallPerf(1,subjects{1})');
    [rho(2,task),p(2,task)] = corr(squeeze(mean(intraSubAvg.CB(subjects{2},:,task),2)), overallPerf(2,subjects{2})');
end

% display results
disp('Pearson correlation between overall performance on recognition tasks and performance on discern pairs')
disp('The first row are the sighted. The second are the blind')
disp('The first column is non-verbal performance. The second is verbal')
disp('Rho: ')
disp(num2str(rho))
disp('P-value')
disp(num2str(p))

figure
for plt = 1:2 % a subplot for nonverbal and verbal
    subplot(1,2,plt)
    hold on
    plot(squeeze(mean(intraSubAvg.S(subjects{1},:,plt),2)),overallPerf(1,subjects{1})','o','MarkerEdgeColor',color(1,:))
    plot(squeeze(mean(intraSubAvg.CB(subjects{2},:,plt),2)),overallPerf(2,subjects{2})','o','MarkerEdgeColor',color(2,:))
    
    hlines=lsline; % weirdness with lsline
    hlines=fliplr(hlines); % need to flip it
    for k = 1:numel(hlines) % then set the colors by hand, argh
        set(hlines(k),'Color',color(k,:))
    end 
    xlabel('Recognition Performance'); ylabel('Discern Pairs Performance')
end


%% End script and turn off diary if writing to text file
if writeToFile
    diary off
end