function Load_All_Results
%copy essential files to each folder
% Get a list of all files and folders in this folder.
files = dir
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% Print folder names to command window.

AllFinalResults={}
t=1
for z=3:size(subFolders, 1)
    cd(subFolders(z).name)
    if exist ('FinalCleanedResults.mat')
        load('FinalCleanedResults.mat','FinalResultsBlast')
        load('FinalCleanedResults.mat','FinalResultsMotif')
        AllFinalResults{t,1}=subFolders(z).name;
        AllFinalResults{t,2}=FinalResultsBlast;
        AllFinalResults{t,3}=FinalResultsMotif;
        t=t+1
    end
    cd ..
end

%Need to import data for locus analysis

%String split for Psuedogene and anotation of full length
%Human

LocusOrder=importdata('C:\Users\GDale\Desktop\Gordon\Human IgHV Gene Order_Listed in order of proximity to D.txt',',')


LocusLocation=LocusOrder.data
LocusOrder=LocusOrder.textdata
LocusOrder=upper(LocusOrder)

% for z=1:length(LocusOrder)
% LocusOrderTemp=strsplit(LocusOrder{z},',')
% LocusOrder{z,1}=LocusOrderTemp{1,1}
% LocusOrder{z,2}=LocusOrderTemp{1,2}
% end


%Check Locus order to see if if 5' or 3' usage and place each count into
%AllFinalResults column 4 and column 5 (5' and 3' respectively)

for x=1:size(AllFinalResults,1)
    FivePrimeCount=0
    ThreePrimeCount=0
    InitialIndex=find(ismember(LocusOrder,upper(AllFinalResults{x,1})))
    AlreadyMatched={}
    for z=1:size(AllFinalResults{x,2},1)
        if isempty(AllFinalResults{x,2}{z,3})==0
            IndexOfHit=find(ismember(LocusOrder,upper(AllFinalResults{x,2}{z,3})))
            AlreadyMatchedTemp=upper(AllFinalResults{x,2}{z,3})
            if exist('IndexOfHit') && ismember(upper(AllFinalResults{x,2}{z,3}),AlreadyMatched)==0
                IndexDiff=IndexOfHit-InitialIndex
                if IndexDiff>0
                    FivePrimeCount=FivePrimeCount+1
                elseif IndexDiff<0
                    ThreePrimeCount=ThreePrimeCount+1
                end
            end
            AlreadyMatched=vertcat(AlreadyMatched,AlreadyMatchedTemp)
            clearvars IndexOfHit
        end
    AllFinalResults{x,4}=FivePrimeCount
    AllFinalResults{x,5}=ThreePrimeCount
    end
end


%Calculate distances between donor Igs and VDJ. Stored in AllFinalResults.
%5' distances stored in column 6, 3' distances stored in column 7
%file. 
for x=1:size(AllFinalResults,1)
    InitialIndex=find(ismember(LocusOrder,upper(AllFinalResults{x,1})))
    AlreadyMatched={}
    k=1
    kk=1
    for z=1:size(AllFinalResults{x,2},1)
        if isempty(AllFinalResults{x,2}{z,3})==0
            IndexOfHit=find(ismember(LocusOrder,upper(AllFinalResults{x,2}{z,3})))
            AlreadyMatchedTemp=upper(AllFinalResults{x,2}{z,3})
            if exist('IndexOfHit') && ismember(upper(AllFinalResults{x,2}{z,3}),AlreadyMatched)==0
                IndexDiff=IndexOfHit-InitialIndex
                if IndexDiff>0
                    AllFinalResults{x,6}(k,1)=LocusLocation(IndexOfHit)-LocusLocation(InitialIndex)
                    k=k+1
                elseif IndexDiff<0
                    AllFinalResults{x,7}(kk,1)=LocusLocation(IndexOfHit)-LocusLocation(InitialIndex)
                    kk=kk+1
                end
            end
            AlreadyMatched=vertcat(AlreadyMatched,AlreadyMatchedTemp)
            clearvars IndexOfHit
        end
    end
end

%Graph Distances of 5' donors (red) and 3' donors (Blue) via histogram
hold on
histogram((unique(vertcat(AllFinalResults{:,6}))),10,'FaceColor','r')
histogram((unique(vertcat(AllFinalResults{:,7}))),10,'FaceColor','b')
hold off
savefig('Cis_Trans.fig')

save('AllDataByIgHV.mat')
end