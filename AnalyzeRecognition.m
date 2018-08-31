%%% Analyze verbal and nonverbal recognition data %%%%

% Script by: Karen Arcos (karcos1@uci.edu) 
% Bedny Lab
% August 2018
clear all; clc;
prompt = 'All data file name?:'; %file name to enter"
f=input(prompt, 's'); %input current file
load(f); %contains structures of all working memory tasks besides recognition (See grabDataMakeStructure.m script used to generate structures)
%span length and performance accuracy column indices from raw recognition data
TrialSpanLengthCol = 7;
lureTypeCol = 11;
isAnswerCorrectCol = 13;
nS = 22; %# sighted participants based on demographics and on those who completed recognition tasks
nCB = 20; %# blind participants based on demographics and on those who completed recognition tasks
n = [nS nCB]; %# participants per group in vector for looping
 tSpans = 4; %total number of spans tested on nonverbal and verbal recognition tasks
triPerSpan = 8; %# of trials per span length
nTriDp = 91; %number of trials tested on nonverbal discern pairs
nTriFwR = tSpans*triPerSpan; %number of trials analyzing for nonverbal and verbal fw recognition
%indices for group, tasks, common spans, and lures
sInd = 1; cbInd = 2; %group indices
nvdp = 4; nvfr = 5; vfr = 6; %task indices
taskIndices = nvdp:vfr; %task indices of each recognition RAW data structure
nvOne = 3; nvTwo = 4; %common span indices for nonverbal recognition spans 5-6
vOne = 1; vTwo = 2; %common span indices for verbal recognition spans 5-6
%nonverbal forward recognition lure indices from REORGANIZED raw data score column (avoids looping extra through lure types)
diffSTwo = 8; diffSThree = 2*diffSTwo; diffSFour = 3*diffSTwo; %values to subtract for span column indices to align to reorganized score column
iChIndNv = [2 4; 9-diffSTwo 16-diffSTwo; 17-diffSThree 22-diffSThree; 27-diffSFour 29-diffSFour]; %identity change lure indices
naIndNv = [1 5 6 7; 11-diffSTwo 12-diffSTwo 14-diffSTwo 15-diffSTwo; 18-diffSThree 19-diffSThree 20-diffSThree 24-diffSThree; 25-diffSFour 26-diffSFour 30-diffSFour 31-diffSFour]; %N/A lure indices
sOOIndNv = [8; 10-diffSTwo; 23-diffSThree; 32-diffSFour]; %SlideOneOver lure indices
sTIndNv = [3;    13-diffSTwo; 21-diffSThree; 28-diffSFour]; %SwapTwo lure indices
%verbal forward recognition lure indices from REORGANIZED raw data score column (avoids looping extra through lure types)
iChIndV = [1 8; 10-diffSTwo 15-diffSTwo;  23-diffSThree 24-diffSThree; 31-diffSFour 32-diffSFour]; %identity change lure indices
naIndV =  [2 3 4 6; 9-diffSTwo 11-diffSTwo 12-diffSTwo 14-diffSTwo; 17-diffSThree 18-diffSThree 20-diffSThree 22-diffSThree; 25-diffSFour 27-diffSFour 29-diffSFour 30-diffSFour]; %N/A lure indices
sOOIndV = [5; 16-diffSTwo; 21-diffSThree; 26-diffSFour]; %SlideOneOver lure indices
sTIndV = [7; 13-diffSTwo; 19-diffSThree; 28-diffSFour]; %SwapTwo lure indices
%lure column indices for lure ACCURACY matrix
aNAInd = 1; aIdChInd = 2; aSOOInd = 3; aSTInd = 4; %lure column indices for lure ACCURACY matrix
%types of subjects we have 
vision = {'S' 'CB'};  %participant labels
data = {'SData' 'CBData'}; %group data labels used for analysis from preloaded file
%initiate structure array to save descriptives per task so 3 structures total
tasks = {nan, nan, nan, 'NonvDP', 'NonvFR', 'VFR'}; %recognition task names; added nan to align indices
groupField = 'SubPercCor'; %percent correct field for each group
recogData = repmat(struct, 1, length(tasks)); %initiate structure to save task raw data and descriptive stats
statsFields = {'averagePercCor', 'stdErr',}; %cell of strings used as descriptive field names based on tasks
delDataS = {'x', 'su11'}; %excluded sighted participants
delDataCB = {'x', 'su3', 'su6', 'su11', 'su12', 'su17'};  %excluded CBsighted participants
   delDataAll = {delDataS delDataCB}; %data to remove (sighted and blind)
for v = 1:length(vision) %loop through blind and sighted groups
functionToEval = strcat(data{v}, '=rmfield(', data{v}, ', delDataAll{v});'); %participant data eval will evaluate
eval(char(functionToEval)); % run the text created in functionToEVal and remove excluded participants
 names = fieldnames(eval(char(data{v}))); %included participant names
 for t = taskIndices %loop through three recognition tasks
    subPerCor = nan(n(v), tSpans); %initiate matrix to save percent correct where participants are rows and spans are columns corresponding to column indices ie col 3 is third span length
aPc = nan(1, tSpans); %initiate matrix to save AVERAGE percents correct across participants where rows are groups and columns are span lengths for one task
sePc = nan(1, length(tSpans)); %initiate matrix to save standard errors of percents correct across participants where columns are span lengths for one task
aPcSe = {'AvePerCor' 'SEPerCor'; repmat(nan,4,1) repmat(nan,4,1)}; %initiating cell to save average percent corrects and standard errors of percents correct per task and group
    spcDiff = nan(n(v), length(v)); %initiate matrix to save percent correct differences where participants are rows and spans 5-6 are columns
    subPerCorLur = nan(tSpans, tSpans, n(v)); %initiate matrix to save percent correct for lures where spans are rows, lure types are columns, and participants are height
hitR = nan(1,n(v)); %initiate row vector to save hit scores averaged across spans
faR = nan(1,n(v)); %initiate row vector to save false alarm scores averaged across spans
dpRow = nan(1,n(v)); %initiate row vector to save d prime scores
for pt = 1:n(v) %loop through participants
functionToEval = strcat(data{v}, '(', num2str(t),').', names{pt}); %participant data eval will evaluate
ptData = eval(char(functionToEval)); % run the text created in functionToEVal and save participant raw data in use
if t == 1
 subPerCor(pt,1) = sum(cell2mat(ptData(:,isAnswerCorrectCol)))/nTriDp; %participant discern pairs percent correct
else %scoring participants for remaining tasks
%adjusting matrix for scoring
[r c] = size(ptData); %size of original sub data
nums = ptData(:,isAnswerCorrectCol); %isolate column for changing class and scoring
	nums = double(cell2mat(nums)); %converts values from logical to double to do math on them
ptData(:,isAnswerCorrectCol) = num2cell(nums); %save values as doubles instead of logical
if r < nTriFwR %when sub not reach max spans tested
ptData=[ptData; repmat({nan},nTriFwR-r,c)]; %add extra rows to cell as needed for max trials tested
useLr = cell2mat(reshape(ptData(:,isAnswerCorrectCol),triPerSpan,tSpans)); %use for lure percents
ptData(r+1:end,isAnswerCorrectCol)=repmat({.5},nTriFwR-r,1); %score correctly on half trials
elseif r >= nTriFwR %when subject exceeds spans tested
ptData = ptData(1:nTriFwR,:); %complete data matrix for scoring by removing rows as needed
r = 0; %use for lure percents
useLr = cell2mat(reshape(ptData(:,isAnswerCorrectCol),triPerSpan,tSpans));
if v==1 & pt ==8
useLr = cell2mat(reshape(ptData(:,isAnswerCorrectCol),triPerSpan,tSpans)) %use for lure percents
end %if v==1
end %if r
usePc=cell2mat(reshape(ptData(:,isAnswerCorrectCol),triPerSpan,tSpans)); %reorganize score column to sum scores per span
subPerCor(pt,:) = sum(usePc)/triPerSpan; %each participant's percent correct per span saved to group matrix
subPerCor(find(subPerCor<.5))=.5; %replace % correctscores lower than  .5 wi/ .5
%lure scores
if t >= nvfr %scores lures for all except nonverbal discern pairs
for s = 1:tSpans %loop going through each span
if t == nvfr %uses nonverbal lure indices
  naInd = naIndNv; iChInd = iChIndNv; sOOInd = sOOIndNv; sTInd = sTIndNv; %use depending on task
elseif t == vfr %uses verbal lure indices
  naInd = naIndV; iChInd = iChIndV; sOOInd = sOOIndV; sTInd = sTIndV; %use depending on task
end %if t
incS = tSpans-((nTriFwR-r)/triPerSpan); %number of incomplete spans
subPerCorLur(s,1:tSpans,pt) = [nansum(useLr(naInd(s,:),s))/(length(naInd)-incS) nansum(useLr(iChInd(s,:),s))/(length(iChInd)-incS) nansum(useLr(sOOInd(s,:),s))/(length(sOOInd)-incS) nansum(useLr(sTInd(s,:),s))/(length(sTInd)-incS)]; %percent correct per lure type and span per participant
end %for s going through each span
%hits and false alarms per participant
hitR(1,pt) = sum(subPerCorLur(:,aNAInd,pt))/(tSpans-incS); %average na lure percent correct scores across spans per participant
if hitR(1,pt)==0 %avoids nan on inversion function
hitR(1,pt)=1/nTriFwR;
elseif hitR(1,pt)==1
hitR(1,pt)=nTriFwR-1/nTriFwR;
end %if hitR
faR(1,pt) = (sum(sum(1-subPerCorLur(:,aIdChInd:aSTInd,pt))))/(tSpans*(tSpans-incS)); %false alarm score found by subtracting total possible correct on each lure (1) and % correct on lures to get total wrong and divide by total lures
if faR(1,pt)==0 %avoids nan on inversion function
faR(1,pt)=1/nTriFwR;
elseif faR(1,pt)==1
faR(1,pt)=(nTriFwR-1)/nTriFwR;
end %if faR
functionToEval = strcat(tasks{t}, '.', vision{v}, '.hits', '= hitR;'); %data eval will evaluate 
% run the text created in functionToEVal
 eval(char(functionToEval)); %save hits per task to structure
functionToEval = strcat(tasks{t}, '.', vision{v}, '.faAl', '= faR;'); %data eval will evaluate 
% run the text created in functionToEVal
 eval(char(functionToEval)); %save false alarms per task to structure

end %if t for lure scoring
functionToEval = strcat(data{v}, '(1,',num2str(t),').', names{pt}, '=[];');
eval(char(functionToEval)); %clear old data from structure
functionToEval = strcat(data{v}, '(',num2str(t),').', names{pt}, '=ptData;'); 
eval(char(functionToEval)); %save new participant data to struct wi eval function
end %if t
end %loop through participants
%descriptive data calculations
%d prime for each task
dpRow(1,:) = norminv(hitR)-norminv(faR); %d prime per participant
%average percent correct across participants
apc = sum(subPerCor)/n(v); %average percent correct per task's span across participants in each group
sePc = std(subPerCor)/sqrt(n(v)); %standard error of percent correct per span across participants in each group
aPcSe{end,1} = apc'; %merge average percent correct & standard error per task so left column is means and right one is standard errors
aPcSe{end,2} = sePc'; %merge average percent correct & standard error per task so left column is means and right one is standard errors
functionToEval = strcat(tasks{t}, '.', vision{v}, '.SubPerCor', '= subPerCor;'); %data eval will evaluate 
% run the text created in functionToEVal
 eval(char(functionToEval)); %saves participants' percents correct per span to tasks structure
functionToEval = strcat(tasks{t}, '.', vision{v}, '.APcSe', '= aPcSe;'); %data eval will evaluate 
% run the text created in functionToEVal
 eval(char(functionToEval)); %save average percents correct and standard error matrix per task to structure
functionToEval = strcat(tasks{t}, '.', vision{v}, '.SubPerCorLur', '= subPerCorLur;'); %data eval will evaluate 
% run the text created in functionToEVal
 eval(char(functionToEval)); %save lure percents correct to structure
functionToEval = strcat(tasks{t}, '.', vision{v}, '.dPr', '= dpRow;'); %data eval will evaluate 
% run the text created in functionToEVal
 eval(char(functionToEval)); %save d primes per task to structure
%one sample t test per group on all recognition tasks
functionToEval = strcat('[', tasks{t}, '.', vision{v}, '.OneSam.H, ', tasks{t}, '.', vision{v}, '.OneSam.P, ', tasks{t}, '.', vision{v}, '.OneSam.CI, ', tasks{t}, '.', vision{v}, '.OneSam.stats] = ttest(', tasks{t}, '.', vision{v}, '.SubPerCor,', num2str(0), ', ''alpha''', ',', num2str(0.01), ');');
eval(char(functionToEval)); %run one sam t test and save to structure
%unpaired samples t test btwn sighted and blind percent correct per span on nv recognition and v recognition
if v == 2 %applies when have data for s and cb on nvfr and vfr
functionToEval = strcat('[', tasks{t}, '.unpSamTPC.H,', tasks{t}, '.unpSamTPC.P,', tasks{t}, '.unpSamTPC.CI,', tasks{t}, '.unpSamTPC.stats] = ttest2(', tasks{t}, '.', vision{v-1}, '.SubPerCor,', tasks{t}, '.', vision{v}, '.SubPerCor);');
eval(char(functionToEval)); %run one sam t test and save to structure
functionToEval = strcat('[', tasks{t}, '.unpSamTPC.H,', tasks{t}, '.unpSamTPC.P,', tasks{t}, '.unpSamTPC.CI,', tasks{t}, '.unpSamTPC.stats] = ttest2(', tasks{t}, '.', vision{v-1}, '.SubPerCor,', tasks{t}, '.', vision{v}, '.SubPerCor);');
eval(char(functionToEval)); %run one sam t test and save to structure
%unpaired samples t test btwn sighted and blind d primes per span on nv recognition and v recognition
functionToEval = strcat('[', tasks{t}, '.unpSamTDPr.H,', tasks{t}, '.unpSamTDPr.P,', tasks{t}, '.unpSamTDPr.CI,', tasks{t}, '.unpSamTDPr.stats] = ttest2(', tasks{t}, '.', vision{v-1}, '.dPr,', tasks{t}, '.', vision{v}, '.dPr);');
eval(char(functionToEval)); %run unpaired sam t test of d pprimes and save to structure
end %if 
end %loop through three recognition tasks
%each participant's difference on verbal and nonverbal percent correct on common spans 5 and 6
functionToEval = strcat('spcDiff(:,1) =', tasks{vfr}, '.', vision{v}, '.SubPerCor(:,', num2str(vOne), ')', '-', tasks{nvfr}, '.', vision{v}, '.SubPerCor(:,', num2str(nvOne), ');'); %subtracting percent corrects on span lengths of 5 for verbal and nonverbal recognition tasks
eval(char(functionToEval)); %save percent correct differences for span length 5 to matrix
functionToEval = strcat('spcDiff(:,2) =', tasks{vfr}, '.', vision{v}, '.SubPerCor(:,', num2str(vTwo), ')', '-', tasks{nvfr}, '.', vision{v}, '.SubPerCor(:,', num2str(nvTwo), ');'); %subtracting percent corrects on span lengths of 6 for verbal and nonverbal recognition tasks
eval(char(functionToEval)); %save percent correct differences for span length 6 to matrix
functionToEval = strcat('stats.spcDiff', vision{v}, '= spcDiff;'); %save each participant's percent correct differences per group
% run the text created in functionToEVal
 eval(char(functionToEval)); %save each group's percent correct differences to structure
functionToEval = strcat('stats.aveSpcDiff', vision{v}, ' = mean(spcDiff);'); %save each group's average percent correct differences
% run the text created in functionToEVal
 eval(char(functionToEval)); %save each group's average percent correct differences to structure
end %loop through blind and sighted groups
%group x task x span (3-way) repeated measures ANOVA btwn groups
taskC = repmat([repmat(nvfr,tSpans,1); repmat(vfr,tSpans,1)],sum(n),1); %task column vector repeating each task index for total spans and participants
groupC = [ones(2*tSpans*nS,1); repmat(cbInd,2*tSpans*nCB,1)]; %column vector repeating each group for total spans and participants
spC=[1:4]'; %span labels (actual span lengths vary from task to task)
spC = repmat(repmat(spC,length(n),1),sum(n),1); %column vector repeating each span label for total tasks and participants
pcCol = nan(tSpans*sum(n),1); %initiating column vector saving all % corrects per span, task, and participant
pcRF = -7; %first percent correct index (updates throughout)
pcRL = 0; %percent correct last index (updates throughout)
for sub = 1:sum(n) %loop going through participants
getPC=nan(2*tSpans,1); %initiating to temporarily saved % corrects
if sub <= nS
 NvPC = strcat('[NonvFR.S.SubPerCor(', num2str(sub), ',:)]'';'); %sighted subject nonverbal %s correct
VPC = strcat('[VFR.S.SubPerCor(', num2str(sub), ',:)]'';'); %sighted subject verbal %s correct
else
NvPC = strcat('[NonvFR.CB.SubPerCor(', num2str(sub-nS), ',:)]'';'); %blind subject nonverbal %s correct
VPC = strcat('[VFR.CB.SubPerCor(', num2str(sub-nS), ',:)]'';'); %blind subject verbal %s correct
end %if sub
getPC(1:tSpans,:) = [eval(char(NvPC))]'; %saves subject %s in use
getPC(tSpans+1:end,:) = [eval(char(VPC))]'; %saves subject %s in use
pcCol(pcRF+2*tSpans:pcRL+2*tSpans,1) = getPC; %save %s to long column
pcRF = pcRF+2*tSpans; pcRL = pcRL+2*tSpans; %new indices for next participant
end %for sub
[p tbl] = anovan(pcCol,{groupC taskC spC},'model','full','varnames',{'groupC','taskC','spC'}, 'display', 'off'); %run anova

%create output file
%save descriptive data structure to output file