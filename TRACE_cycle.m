function TRACE_cycle(MotifLength, NumMuts, InputFileHandle, BlastDB, AnalysisName,Path_BlastDB)
names = dir(InputFileHandle)
names = {names.name}
tic

for Cycle=1:max(size(names))
    Cycle
    inputFasta=names{Cycle};
    TRACE_noGaps_1000(inputFasta, MotifLength, NumMuts, BlastDB, AnalysisName,Path_BlastDB)
    
    load((sprintf('%s_TRACE_%s.mat', AnalysisName, inputFasta)), 'sampleZScore');
  
    Dirname=(sprintf('%s', names{Cycle}));
    Dirname=erase(Dirname, '.fas');
    mkdir(Dirname);  
    movefile('*TRACE*.mat', Dirname);
    
    if exist('sampleZScore')>0
    movefile('*.fasta',Dirname);
    else end
    
    clearvars -except MotifLength NumMuts InputFileHandle names BlastDB AnalysisName
end

HomeFolder=pwd

% Get a list of all files and folders in this folder.
files = dir
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% Print folder names to command window.

for k = 1 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end
subFolders={subFolders.name};
subFolders=subFolders(1,3:length(subFolders));

clearvars -except subFolders AnalysisName HomeFolder
for Cycle=1:max(size(subFolders))
cd(sprintf('%s', subFolders{Cycle}))
if exist((sprintf('%s_TRACE_%s%s.mat', AnalysisName,subFolders{Cycle},'.fas')),'file')==2
load((sprintf('%s_TRACE_%s%s.mat', AnalysisName,subFolders{Cycle},'.fas')), 'sampleZScore');
load((sprintf('%s_TRACE_%s%s.mat', AnalysisName,subFolders{Cycle},'.fas')), 'stoufferZ');
load((sprintf('%s_TRACE_%s%s.mat', AnalysisName,subFolders{Cycle},'.fas')), 'originalstoufferZ');
load((sprintf('%s_TRACE_%s%s.mat', AnalysisName,subFolders{Cycle},'.fas')), 'originalsampleZScore');

if exist('sampleZScore')>0
sampleZScore=transpose(sampleZScore);
FinalsampleZScoreCombined{:,Cycle}=sampleZScore;
stoufferZCombined{Cycle,1}=stoufferZ;
originalsampleZScore=transpose(originalsampleZScore)
FinaloriginalsampleZScoreCombined{:,Cycle}=originalsampleZScore;
originalstoufferZCombined{Cycle,1}=originalstoufferZ;
end

else
end

cd(HomeFolder)

clearvars -except subFolders FinalsampleZScoreCombined FinalsampleZScore_decompressedCombined stoufferZCombined stoufferZ_decompressedCombined AnalysisName HomeFolder originalstoufferZCombined FinaloriginalsampleZScoreCombined
end
%------------------------Phase2 Calc----------------------------------
if exist('FinaloriginalsampleZScoreCombined')==1
Final_stoufferZCombined=vertcat(stoufferZCombined{:});
for t=1:max(size(Final_stoufferZCombined))
if isempty(Final_stoufferZCombined(t))==0
    if isfinite(Final_stoufferZCombined(t))==1
        Final_stoufferZCombined2{t}=Final_stoufferZCombined(t)
    elseif isfinite(Final_stoufferZCombined(t))==0
        Final_stoufferZCombined2{t}=[]
    end
else
end
end
Final_stoufferZCombined=vertcat(Final_stoufferZCombined2{:})
ZScoreCount=vertcat(FinalsampleZScoreCombined{:});
ZScoreCount=ZScoreCount(isfinite(ZScoreCount));
ZScoreCount=max(size(ZScoreCount));
z=0;

for Cycle2=1:max(size(FinalsampleZScoreCombined(1,:)));
z=z+1;
if isempty(FinalsampleZScoreCombined{1,Cycle2})<1;
FinalsampleZScoreCombined2=cell2mat(FinalsampleZScoreCombined(1,Cycle2));
FinalsampleZScoreCombined2=FinalsampleZScoreCombined2(isfinite(FinalsampleZScoreCombined2));
Weights(z,1)=(max(size(FinalsampleZScoreCombined2)))/(ZScoreCount);
WeightsSquared(z,1)=(Weights(z,1))*(Weights(z,1));
if isfinite(stoufferZCombined{z,1})==0;
stoufferZCombined{z,1}=sum(FinalsampleZScoreCombined2)/((max(size(FinalsampleZScoreCombined2)))^0.5);
else end
WeightedStouffers(z,1)=(stoufferZCombined{z,1})*(Weights(z,1));
clearvars FinalsampleZScoreCombined2

elseif isempty(FinalsampleZScoreCombined{1,Cycle2})==1 ;
end
end

StouffersZTrend_Phase2=(sum(WeightedStouffers))/(sqrt(sum(WeightsSquared)))


%Remove Inf from analysis 
Final_originalstoufferZCombined=vertcat(originalstoufferZCombined{:});
for t=1:max(size(Final_originalstoufferZCombined))
if isempty(Final_originalstoufferZCombined(t))==0
    if isfinite(Final_originalstoufferZCombined(t))==1
        Final_originalstoufferZCombined2{t}=Final_originalstoufferZCombined(t)
    elseif isfinite(Final_originalstoufferZCombined(t))==0
        Final_originalstoufferZCombined2{t}=[]
    end
else
end
end
Final_originalstoufferZCombined=vertcat(Final_originalstoufferZCombined2{:})

%------------------------Phase1 Calc----------------------------------
originalZScoreCount=vertcat(FinaloriginalsampleZScoreCombined{:});
originalZScoreCount=originalZScoreCount(isfinite(originalZScoreCount));
originalZScoreCount=max(size(originalZScoreCount));
z=0;

for Cycle2=1:max(size(FinaloriginalsampleZScoreCombined(1,:)));
z=z+1;
if isempty(FinaloriginalsampleZScoreCombined{1,Cycle2})<1;
FinaloriginalsampleZScoreCombined2=cell2mat(FinaloriginalsampleZScoreCombined(1,Cycle2));
FinaloriginalsampleZScoreCombined2=FinaloriginalsampleZScoreCombined2(isfinite(FinaloriginalsampleZScoreCombined2));
originalWeights(z,1)=(max(size(FinaloriginalsampleZScoreCombined2)))/(originalZScoreCount);
originalWeightsSquared(z,1)=(originalWeights(z,1))*(originalWeights(z,1));
if isfinite(originalstoufferZCombined{z,1})==0;
originalstoufferZCombined{z,1}=sum(FinaloriginalsampleZScoreCombined2)/((max(size(FinaloriginalsampleZScoreCombined2)))^0.5);
else end
originalWeightedStouffers(z,1)=(originalstoufferZCombined{z,1})*(originalWeights(z,1));
clearvars FinaloriginalsampleZScoreCombined2

elseif isempty(FinaloriginalsampleZScoreCombined{1,Cycle2})==1 ;
end
end

StouffersZTrend_Phase1=(sum(originalWeightedStouffers))/(sqrt(sum(originalWeightsSquared)))
runtime=toc
save(sprintf('TRACE_%s.mat', AnalysisName))
end
end
