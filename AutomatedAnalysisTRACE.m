function AutomatedAnalysisTRACE(InputFileHandle,InRange, NumMuts, NumReps, BlastDB, AnalysisName, species, MaxNumMuts,Path_RequiredFiles, Path_BlastDB)
% This function serves as the main function for executing TRACE. 
%% Input Arguments:
% InputFileHandle - File extension for fasta input. Usually '*.fas',
% '*.fa', or '*.fasta'. Asterisk required by MATLAB and serves as a
% wildcard. 

% InRange - the Motif size to be used in the TRACE pipeline. As configured
% in the paper this value is 38.

% NumMuts - the number of mutations within a given InRange value to qualify
% as a mutation cluster. As configured in the paper this value is 8.

% NumReps - the number of background datasets to generate. As configured in
% the paper this value is 10.

% BlastDB - the FASTA file for which a local BLAST database is generated
% from and which will be used for the subsequent TRACE analysis. As
% configured in the paper this was the hg38 Human Reference Genome. 

% AnalysisName - User defined name for the TRACE run. Will be stored in
% output .mat file

% species - can be either 'Human' or 'Mouse'. Controls chromosome number 
% and gene coordinates in the TRACE_PostProcessing step

% MaxNumMuts - controls the number of bins in the output Histogram figures
% created in the Mutation analysis steps. As configred in the paper, this
% value is 30.

% Path_RequiredFiles - Path to the RequiredFiles Folder. *Must not include
% a space in path name*

% Path_BlastDB - Path to the BlastDB folder where the BlastDB fasta file is
% located *Must not include a space in path name*

%% Overview:
% First background data sets are generated, and the test sequence is placed 
% in the 'Original Seqs' folder. Then necessary TRACE files are copied into
% the created directories. The script then cycles through each folder
% executing the TRACE pipeline. At the end, the data is loaded together,
% and the pipeline analyzes the TRACE hits, extracts the nucleotide
% sequence of the search space (as configured in our paper this is the
% human genome), and quantitates the number of mutations explained by the
% TRACE identified template. The background data and original data are
% placed in bins, frequencies are calculated, and the final data is
% extracted. Final data is placed into the variable FinalResults. This
% includes data on the BLAST hit and on the Motif used to generate this
% blast hit.

%% Output Format
% Data in FinalResults is structured as follows:
% Col 1 -> raw BLAST output. Presented in standard BLAST ouput. Sequence
% Header is represented by nucleotide sequence (presented as numbers (A->1,
% C->2, G->3, T->4)) followed by an "_", the sequence number, an "_", and
% starting position in the sequence. The final number is an identifier.
% Col 2-> Exon/Intron/Intergenic region annotation
% Col 3-> Gene name (if applicable)
% Col 4-> Strand orientation
% Col 5-> Chromosome 
% Col 6-> Length of BLAST hit
% Col 7-> Nucleotide sequence, translated version on that in Col 1.
% Col 8-> 0's and 1's representing sites of mutations in the motif. 1's
% denote mutations, 0's are germline. 
% Col 9-> Number of mutations in the motif
% Col 10-> BLAST percent identity with output
% Col 11-> BLAST length
% Col 12-> Number of mutations explained by the BLAST hit

%% Data formating and Background data set creation
ModelingMutationPos_CoordinateScramble(InputFileHandle, InRange, NumMuts, NumReps)

mkdir('Original Seqs')
movefile(InputFileHandle,'Original Seqs')

%copy essential files to each folder
% Get a list of all files and folders in this folder.
files = dir
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% Print folder names to command window.
for z=3:size(subFolders, 1)
    copyfile('fastawriteBetter.m',subFolders(z).name)
    copyfile('TRACE_cycle.m',subFolders(z).name)
    copyfile('TRACE_noGaps_1000.m',subFolders(z).name)
    copyfile('TRACE_PostProcessing.m',subFolders(z).name)
end

%% TRACE analysis of each folder
for z=3:size(subFolders, 1)
    disp(subFolders(z).name)
    cd(subFolders(z).name)
    if exist('MAT Files', 'file')==7
        rmdir('MAT Files', 's')
    end
    TRACE_cycle(InRange,NumMuts,InputFileHandle,BlastDB,AnalysisName,Path_BlastDB);
    TRACE_PostProcessing(species, Path_RequiredFiles)
    cd ..
end

%% Mutation Analysis
MutationCoverageAnalysis_PostTRACEv2('Background',InputFileHandle,MaxNumMuts);
MutationCoverageAnalysis_PostTRACEv2('Original',InputFileHandle,MaxNumMuts);

%% Background Mutation Analysis
load('ConcatenatedResults_Background.mat')

clearvars -except TotalCompiledStorage PercentMatchList polymorphismPos A ConcatenatedUniquePercentMatchingMotifs BlastDB MaxNumMuts Reference Path_BlastDB
for H=1:size(PercentMatchList,2)
    try
        K=1;
        for y=1:size(TotalCompiledStorage,1)
            temp1=strsplit(TotalCompiledStorage{y,1},',');
            if str2double(temp1{3})==str2double(PercentMatchList{H})
                temp2=strsplit(temp1{1,1},'_');
                length(temp2{1,1});
                qstart=str2num(temp1{7});
                qend=str2num(temp1{8});
                z=1;
                for x=qstart:qend
                    temp3(z)=str2num(temp2{1,1}(x));
                    z=z+1;
                end
                PercentMatchingMotifs{K,1}=int2nt(temp3);
                ResortedCompiledStorage(K,:)=TotalCompiledStorage(y,:);
                motifmatchsizes{K,1}=temp1{4};
                clearvars('temp3','temp2','temp1')
                K=K+1;
            end
        end
        
        UniquePercentMatchingMotifs=PercentMatchingMotifs;
        Uniquemotifmatchsizes=motifmatchsizes;
        
        %Get Positions of Matching Motifs in Fasta File
        K=1;
        for Z=1:size(UniquePercentMatchingMotifs(:,:),1)
            M=0;
            V=1;
            while M==0
                MotifMatchTemp=strfind(A(V).Sequence, PercentMatchingMotifs(Z));
                if isempty(MotifMatchTemp)==0
                    MotifMatchCoordinates{K,1}=V;
                    MotifMatchCoordinates{K,2}=MotifMatchTemp;
                    MotifMatchCoordinates{K,3}=MotifMatchTemp+(length(UniquePercentMatchingMotifs{Z,1})-1);
                    K=K+1;
                    M=1;
                end
                V=V+1;
            end
            clearvars('M','V','z');
        end
        
        %Extract string of PolymorphismPos to compare next to PerfectMotifMatch
        for c=1:size(MotifMatchCoordinates(:,:),1)
            y=1;
            for x=MotifMatchCoordinates{c,2}:MotifMatchCoordinates{c,3};
                MotifPolymorphismTemp(1,y)=polymorphismPos(MotifMatchCoordinates{c,1},x);
                y=y+1;
            end
            MotifPolymorphismTemp=double(MotifPolymorphismTemp);
            MotifPolymorphismTemp2=string(MotifPolymorphismTemp);
            MotifPolymorphismTemp2=strjoin(MotifPolymorphismTemp2);
            MotifPolymorphismTemp2= MotifPolymorphismTemp2{1}(find(~isspace(MotifPolymorphismTemp2)));
            UniquePercentMatchingMotifs{c,2}=MotifPolymorphismTemp2;
            UniquePercentMatchingMotifs{c,3}=sum(MotifPolymorphismTemp);
            UniquePercentMatchingMotifs{c,4}=PercentMatchList{H};
            clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
        end
        
        if H>1;
            UniquePercentMatchingMotifs2=vertcat(UniquePercentMatchingMotifs2, UniquePercentMatchingMotifs);
            Uniquemotifmatchsizes2=vertcat(Uniquemotifmatchsizes2, Uniquemotifmatchsizes);
            ResortedCompiledStorage2=vertcat(ResortedCompiledStorage2, ResortedCompiledStorage);
        else
            UniquePercentMatchingMotifs2=UniquePercentMatchingMotifs;
            Uniquemotifmatchsizes2=Uniquemotifmatchsizes;
            ResortedCompiledStorage2=ResortedCompiledStorage;
        end
        clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts UniquePercentMatchingMotifs2 Uniquemotifmatchsizes2 ResortedCompiledStorage2 BlastDB ConcatenatedUniquePercentMatchingMotifs MaxNumMuts Reference Path_BlastDB
    catch
        try
            H=H+1
            K=1;
            for y=1:size(TotalCompiledStorage,1)
                temp1=strsplit(TotalCompiledStorage{y,1},',');
                if str2double(temp1{3})==str2double(PercentMatchList{H})
                    temp2=strsplit(temp1{1,1},'_');
                    length(temp2{1,1});
                    qstart=str2num(temp1{7});
                    qend=str2num(temp1{8});
                    z=1;
                    for x=qstart:qend
                        temp3(z)=str2num(temp2{1,1}(x));
                        z=z+1;
                    end
                    PercentMatchingMotifs{K,1}=int2nt(temp3);
                    ResortedCompiledStorage(K,:)=TotalCompiledStorage(y,:);
                    motifmatchsizes{K,1}=temp1{4};
                    clearvars('temp3','temp2','temp1')
                    K=K+1;
                end
            end
            
            UniquePercentMatchingMotifs=PercentMatchingMotifs;
            Uniquemotifmatchsizes=motifmatchsizes;
            
            %Get Positions of Matching Motifs in Fasta File
            K=1;
            for Z=1:size(UniquePercentMatchingMotifs(:,:),1)
                M=0;
                V=1;
                while M==0
                    MotifMatchTemp=strfind(A(V).Sequence, PercentMatchingMotifs(Z));
                    if isempty(MotifMatchTemp)==0
                        MotifMatchCoordinates{K,1}=V;
                        MotifMatchCoordinates{K,2}=MotifMatchTemp;
                        MotifMatchCoordinates{K,3}=MotifMatchTemp+(length(UniquePercentMatchingMotifs{Z,1})-1);
                        K=K+1;
                        M=1;
                    end
                    V=V+1;
                end
                clearvars('M','V','z');
            end
            
            %Extract string of PolymorphismPos to compare next to PerfectMotifMatch
            for c=1:size(MotifMatchCoordinates(:,:),1)
                y=1;
                for x=MotifMatchCoordinates{c,2}:MotifMatchCoordinates{c,3};
                    MotifPolymorphismTemp(1,y)=polymorphismPos(MotifMatchCoordinates{c,1},x);
                    y=y+1;
                end
                MotifPolymorphismTemp=double(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=string(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=strjoin(MotifPolymorphismTemp2);
                MotifPolymorphismTemp2= MotifPolymorphismTemp2{1}(find(~isspace(MotifPolymorphismTemp2)));
                UniquePercentMatchingMotifs{c,2}=MotifPolymorphismTemp2;
                UniquePercentMatchingMotifs{c,3}=sum(MotifPolymorphismTemp);
                UniquePercentMatchingMotifs{c,4}=PercentMatchList{H};
                clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
            end
            
            if H>1;
                UniquePercentMatchingMotifs2=vertcat(UniquePercentMatchingMotifs2, UniquePercentMatchingMotifs);
                Uniquemotifmatchsizes2=vertcat(Uniquemotifmatchsizes2, Uniquemotifmatchsizes);
                ResortedCompiledStorage2=vertcat(ResortedCompiledStorage2, ResortedCompiledStorage);
            else
                UniquePercentMatchingMotifs2=UniquePercentMatchingMotifs;
                Uniquemotifmatchsizes2=Uniquemotifmatchsizes;
                ResortedCompiledStorage2=ResortedCompiledStorage;
            end
            clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts UniquePercentMatchingMotifs2 Uniquemotifmatchsizes2 ResortedCompiledStorage2 BlastDB ConcatenatedUniquePercentMatchingMotifs MaxNumMuts Reference Path_BlastDB
        catch
            H=H+2
            K=1;
            for y=1:size(TotalCompiledStorage,1)
                temp1=strsplit(TotalCompiledStorage{y,1},',');
                if str2double(temp1{3})==str2double(PercentMatchList{H})
                    temp2=strsplit(temp1{1,1},'_');
                    length(temp2{1,1});
                    qstart=str2num(temp1{7});
                    qend=str2num(temp1{8});
                    z=1;
                    for x=qstart:qend
                        temp3(z)=str2num(temp2{1,1}(x));
                        z=z+1;
                    end
                    PercentMatchingMotifs{K,1}=int2nt(temp3);
                    ResortedCompiledStorage(K,:)=TotalCompiledStorage(y,:);
                    motifmatchsizes{K,1}=temp1{4};
                    clearvars('temp3','temp2','temp1')
                    K=K+1;
                end
            end
            
            UniquePercentMatchingMotifs=PercentMatchingMotifs;
            Uniquemotifmatchsizes=motifmatchsizes;
            
            %Get Positions of Matching Motifs in Fasta File
            K=1;
            for Z=1:size(UniquePercentMatchingMotifs(:,:),1)
                M=0;
                V=1;
                while M==0
                    MotifMatchTemp=strfind(A(V).Sequence, PercentMatchingMotifs(Z));
                    if isempty(MotifMatchTemp)==0
                        MotifMatchCoordinates{K,1}=V;
                        MotifMatchCoordinates{K,2}=MotifMatchTemp;
                        MotifMatchCoordinates{K,3}=MotifMatchTemp+(length(UniquePercentMatchingMotifs{Z,1})-1);
                        K=K+1;
                        M=1;
                    end
                    V=V+1;
                end
                clearvars('M','V','z');
            end
            
            %Extract string of PolymorphismPos to compare next to PerfectMotifMatch
            for c=1:size(MotifMatchCoordinates(:,:),1)
                y=1;
                for x=MotifMatchCoordinates{c,2}:MotifMatchCoordinates{c,3};
                    MotifPolymorphismTemp(1,y)=polymorphismPos(MotifMatchCoordinates{c,1},x);
                    y=y+1;
                end
                MotifPolymorphismTemp=double(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=string(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=strjoin(MotifPolymorphismTemp2);
                MotifPolymorphismTemp2= MotifPolymorphismTemp2{1}(find(~isspace(MotifPolymorphismTemp2)));
                UniquePercentMatchingMotifs{c,2}=MotifPolymorphismTemp2;
                UniquePercentMatchingMotifs{c,3}=sum(MotifPolymorphismTemp);
                UniquePercentMatchingMotifs{c,4}=PercentMatchList{H};
                clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
            end
            
            if H>1;
                UniquePercentMatchingMotifs2=vertcat(UniquePercentMatchingMotifs2, UniquePercentMatchingMotifs);
                Uniquemotifmatchsizes2=vertcat(Uniquemotifmatchsizes2, Uniquemotifmatchsizes);
                ResortedCompiledStorage2=vertcat(ResortedCompiledStorage2, ResortedCompiledStorage);
            else
                UniquePercentMatchingMotifs2=UniquePercentMatchingMotifs;
                Uniquemotifmatchsizes2=Uniquemotifmatchsizes;
                ResortedCompiledStorage2=ResortedCompiledStorage;
            end
            clearvars -except Path_BlastDB PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts UniquePercentMatchingMotifs2 Uniquemotifmatchsizes2 ResortedCompiledStorage2 BlastDB ConcatenatedUniquePercentMatchingMotifs MaxNumMuts Reference
        end
    end
end
UniquePercentMatchingMotifs2(:,5)=Uniquemotifmatchsizes2;

%Define Number of Mutations Explained by TRACE results
MotifsUsedByTRACE(:,1)=UniquePercentMatchingMotifs2(:,2);
MotifsUsedByTRACE(:,2)=UniquePercentMatchingMotifs2(:,1)
for x=1:size(MotifsUsedByTRACE,1)
    Temp={};
    for ggg=1:size(MotifsUsedByTRACE{x,1},2)
        Temp2=MotifsUsedByTRACE{x,1}(ggg);
        Temp=horzcat(Temp,Temp2);
    end
    Temp=transpose(str2num(str2mat(Temp)));
    for y=1:length(Temp)
        if Temp(y)==0
            Temp(y)=2;
        end
        if Temp(y)==1
            Temp(y)=3;
        end
    end
    TranslatedMotifsUsedByTRACE{x,1}=int2nt(Temp);
    TranslatedMotifsUsedByTRACE_Header{x,1}=sprintf('Motif_%d',x)
end

for x=1:size(MotifsUsedByTRACE,1)
    TranslatedMotifsUsedByTRACE0{x,1}=TranslatedMotifsUsedByTRACE{x,1};
    TranslatedMotifsUsedByTRACESeq{x,1}=MotifsUsedByTRACE{x,2};
end


fastawrite('MotifsUsedByTRACE_MutationLocation_Background.fas',TranslatedMotifsUsedByTRACE_Header,TranslatedMotifsUsedByTRACE0)
fastawrite('MotifsUsedByTRACE_Sequences_Background.fas',TranslatedMotifsUsedByTRACE_Header,TranslatedMotifsUsedByTRACESeq)
save('MotifsUsedByTRACE_Analysis_Background.mat')

A=fastaread('MotifsUsedByTRACE_Sequences_Background.fas');
Reference=fastaread(sprintf('%s%s','C:\Users\GDale\Desktop\Gordon\TRACE_Databases\Files\',BlastDB));
ReferenceRegions={}
for x=1:size(A,1)
    x
    fastawrite('distIndex.fasta',A(x))
    [~,blastOutput]=system(sprintf('blastn.exe -db "%s%s" -query distIndex.fasta -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -word_size 11 -max_hsps 1 -max_target_seqs 1 -num_threads 2', Path_BlastDB, BlastDB));
    BlastTemp1=strsplit(blastOutput,',')
    delete distIndex.fasta
    pause(0.25)
    if strcmp(strtrim(BlastTemp1{13}),'plus');
        CoordinateStart=str2double(BlastTemp1{9});
        CoordinateEnd=str2num(BlastTemp1{10});
        for y=1:size(Reference,1)
            if strcmp(Reference(y).Header,BlastTemp1{2})==1
                k=1;
                for z=CoordinateStart:CoordinateEnd
                    ReferenceRegionsTemp(k)=Reference(y).Sequence(z);
                    k=k+1;
                end
                ReferenceRegionsTemp2{1}=(ReferenceRegionsTemp);
                ReferenceRegions{x,1}=ReferenceRegionsTemp2{1};
                clearvars ReferenceRegionsTemp ReferenceRegionsTemp2
            end
        end
    elseif strcmp(strtrim(BlastTemp1{13}),'minus');
        CoordinateStart=str2double(BlastTemp1{10});
        CoordinateEnd=str2num(BlastTemp1{9});
        for y=1:size(Reference,1)
            if strcmp(Reference(y).Header,BlastTemp1{2})==1
                k=1;
                for z=CoordinateStart:CoordinateEnd
                    ReferenceRegionsTemp(k)=Reference(y).Sequence(z);
                    k=k+1;
                end
                ReferenceRegionsTemp2{1}=seqrcomplement(ReferenceRegionsTemp);
                ReferenceRegions{x,1}=ReferenceRegionsTemp2{1};;
                clearvars ReferenceRegionsTemp ReferenceRegionsTemp2
            end
        end
    end
    clearvars BlastTemp1
end
ReferenceRegions=upper(ReferenceRegions); %control for sequences that return lowercase
load('MotifsUsedByTRACE_Analysis_Background.mat','MotifsUsedByTRACE')

MutationExplainedBlastMap={}
for x=1:size(MotifsUsedByTRACE)
    if length((ReferenceRegions{x}))-length(MotifsUsedByTRACE{x,2})~=0
        [Score,Alignment]=nwalign(ReferenceRegions{x},MotifsUsedByTRACE{x,2})
        ReferenceRegions{x}=Alignment(1,:);
        MotifsUsedByTRACE{x,2}=Alignment(3,:);
        if length(str2num(str2mat(MotifsUsedByTRACE{x,1})))-length(MotifsUsedByTRACE{x,2})~=0
            Intermediate={};
            for ggg=1:size(MotifsUsedByTRACE{x,1},2)
                Intermediate2=MotifsUsedByTRACE{x,1}(ggg);
                Intermediate=horzcat(Intermediate,Intermediate2);
            end
            Intermediate=transpose(str2num(str2mat(Intermediate)));
        end
    end
    testDataset=vertcat(ReferenceRegions{x},MotifsUsedByTRACE{x,2});
    polymorphismPos=testDataset~=repmat(testDataset(1,:), size(testDataset, 1), 1);
    polymorphismPos(testDataset=='-')=0;
    if strfind(ReferenceRegions{x},'-')~=0
        GapsToCorrect=strfind(ReferenceRegions{x},'-')
        for zz=1:length(GapsToCorrect)
            polymorphismPos(2,GapsToCorrect(zz))=0
        end
    end
    if strfind(MotifsUsedByTRACE{x,2},'-')~=0
        GapsToCorrect=strfind(MotifsUsedByTRACE{x,2},'-')
        for zz=1:length(GapsToCorrect)
            polymorphismPos(2,GapsToCorrect(zz))=0
            Intermediate=[Intermediate(:,1:GapsToCorrect(zz)-1) 0 Intermediate(:,GapsToCorrect(zz):end)];
            MotifsUsedByTRACE{x,1}=mat2str(Intermediate)
        end
    end
    
    %hotfix
    FixedMotif={};
    for ggg=1:size(MotifsUsedByTRACE{x,1},2)
        Intermediate2=MotifsUsedByTRACE{x,1}(ggg);
        FixedMotif=horzcat(FixedMotif,Intermediate2);
    end
    FixedMotif=transpose(str2num(str2mat(FixedMotif)));
    
    MutationExplainedBlastMap{x,1}=polymorphismPos(2,:)-FixedMotif;
    MutationExplainedPerHit{x,1}=sum(MutationExplainedBlastMap{x,1}(:) == -1);
    clearvars testDataset polymorphismPos Intermediate
end
UniquePercentMatchingMotifs2(:,6)=MutationExplainedPerHit

%Graph 2DHeatmap for Explained Mutations
for x=1:size(UniquePercentMatchingMotifs2,1)
    ExpNumMuts(x)=UniquePercentMatchingMotifs2{x,6};
    PercentMatch(x)=str2num(UniquePercentMatchingMotifs2{x,4});
end

X=[ExpNumMuts;PercentMatch];
X=transpose(X);
hist3(X,'Ctrs',{0:1:MaxNumMuts 80:2:100},'CdataMode','auto')
xlabel('NumMuts')
ylabel('PercentMatch')
colorbar
view(2)

savefig(sprintf('%s_%s','HistogramHitsCorrected','Background'))

save('Backgrounddata_by_matchlength.mat')
clearvars -except BlastDB MaxNumMuts Reference Path_BlastDB

%% Original Mutation Analysis
load('ConcatenatedResults_Original.mat')

clearvars -except TotalCompiledStorage PercentMatchList polymorphismPos A ConcatenatedUniquePercentMatchingMotifs BlastDB MaxNumMuts Reference Path_BlastDB
for H=1:size(PercentMatchList,2)
    try
        K=1;
        for y=1:size(TotalCompiledStorage,1)
            temp1=strsplit(TotalCompiledStorage{y,1},',');
            if str2double(temp1{3})==str2double(PercentMatchList{H})
                temp2=strsplit(temp1{1,1},'_');
                length(temp2{1,1});
                qstart=str2num(temp1{7});
                qend=str2num(temp1{8});
                z=1;
                for x=qstart:qend
                    temp3(z)=str2num(temp2{1,1}(x));
                    z=z+1;
                end
                PercentMatchingMotifs{K,1}=int2nt(temp3);
                ResortedCompiledStorage(K,:)=TotalCompiledStorage(y,:);
                motifmatchsizes{K,1}=temp1{4};
                clearvars('temp3','temp2','temp1')
                K=K+1;
            end
        end
        
        UniquePercentMatchingMotifs=PercentMatchingMotifs;
        Uniquemotifmatchsizes=motifmatchsizes;
        
        %Get Positions of Matching Motifs in Fasta File
        K=1;
        for Z=1:size(UniquePercentMatchingMotifs(:,:),1)
            M=0;
            V=1;
            while M==0
                MotifMatchTemp=strfind(A(V).Sequence, PercentMatchingMotifs(Z));
                if isempty(MotifMatchTemp)==0
                    MotifMatchCoordinates{K,1}=V;
                    MotifMatchCoordinates{K,2}=MotifMatchTemp;
                    MotifMatchCoordinates{K,3}=MotifMatchTemp+(length(UniquePercentMatchingMotifs{Z,1})-1);
                    K=K+1;
                    M=1;
                end
                V=V+1;
            end
            clearvars('M','V','z');
        end
        
        %Extract string of PolymorphismPos to compare next to PerfectMotifMatch
        for c=1:size(MotifMatchCoordinates(:,:),1)
            y=1;
            for x=MotifMatchCoordinates{c,2}:MotifMatchCoordinates{c,3};
                MotifPolymorphismTemp(1,y)=polymorphismPos(MotifMatchCoordinates{c,1},x);
                y=y+1;
            end
            MotifPolymorphismTemp=double(MotifPolymorphismTemp);
            MotifPolymorphismTemp2=string(MotifPolymorphismTemp);
            MotifPolymorphismTemp2=strjoin(MotifPolymorphismTemp2);
            MotifPolymorphismTemp2= MotifPolymorphismTemp2{1}(find(~isspace(MotifPolymorphismTemp2)));
            UniquePercentMatchingMotifs{c,2}=MotifPolymorphismTemp2;
            UniquePercentMatchingMotifs{c,3}=sum(MotifPolymorphismTemp);
            UniquePercentMatchingMotifs{c,4}=PercentMatchList{H};
            clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
        end
        
        if H>1;
            UniquePercentMatchingMotifs2=vertcat(UniquePercentMatchingMotifs2, UniquePercentMatchingMotifs);
            Uniquemotifmatchsizes2=vertcat(Uniquemotifmatchsizes2, Uniquemotifmatchsizes);
            ResortedCompiledStorage2=vertcat(ResortedCompiledStorage2, ResortedCompiledStorage);
        else
            UniquePercentMatchingMotifs2=UniquePercentMatchingMotifs;
            Uniquemotifmatchsizes2=Uniquemotifmatchsizes;
            ResortedCompiledStorage2=ResortedCompiledStorage;
        end
        clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts UniquePercentMatchingMotifs2 Uniquemotifmatchsizes2 ResortedCompiledStorage2 BlastDB ConcatenatedUniquePercentMatchingMotifs MaxNumMuts Reference Path_BlastDB
    catch
        try
            H=H+1
            K=1;
            for y=1:size(TotalCompiledStorage,1)
                temp1=strsplit(TotalCompiledStorage{y,1},',');
                if str2double(temp1{3})==str2double(PercentMatchList{H})
                    temp2=strsplit(temp1{1,1},'_');
                    length(temp2{1,1});
                    qstart=str2num(temp1{7});
                    qend=str2num(temp1{8});
                    z=1;
                    for x=qstart:qend
                        temp3(z)=str2num(temp2{1,1}(x));
                        z=z+1;
                    end
                    PercentMatchingMotifs{K,1}=int2nt(temp3);
                    ResortedCompiledStorage(K,:)=TotalCompiledStorage(y,:);
                    motifmatchsizes{K,1}=temp1{4};
                    clearvars('temp3','temp2','temp1')
                    K=K+1;
                end
            end
            
            UniquePercentMatchingMotifs=PercentMatchingMotifs;
            Uniquemotifmatchsizes=motifmatchsizes;
            
            %Get Positions of Matching Motifs in Fasta File
            K=1;
            for Z=1:size(UniquePercentMatchingMotifs(:,:),1)
                M=0;
                V=1;
                while M==0
                    MotifMatchTemp=strfind(A(V).Sequence, PercentMatchingMotifs(Z));
                    if isempty(MotifMatchTemp)==0
                        MotifMatchCoordinates{K,1}=V;
                        MotifMatchCoordinates{K,2}=MotifMatchTemp;
                        MotifMatchCoordinates{K,3}=MotifMatchTemp+(length(UniquePercentMatchingMotifs{Z,1})-1);
                        K=K+1;
                        M=1;
                    end
                    V=V+1;
                end
                clearvars('M','V','z');
            end
            
            %Extract string of PolymorphismPos to compare next to PerfectMotifMatch
            for c=1:size(MotifMatchCoordinates(:,:),1)
                y=1;
                for x=MotifMatchCoordinates{c,2}:MotifMatchCoordinates{c,3};
                    MotifPolymorphismTemp(1,y)=polymorphismPos(MotifMatchCoordinates{c,1},x);
                    y=y+1;
                end
                MotifPolymorphismTemp=double(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=string(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=strjoin(MotifPolymorphismTemp2);
                MotifPolymorphismTemp2= MotifPolymorphismTemp2{1}(find(~isspace(MotifPolymorphismTemp2)));
                UniquePercentMatchingMotifs{c,2}=MotifPolymorphismTemp2;
                UniquePercentMatchingMotifs{c,3}=sum(MotifPolymorphismTemp);
                UniquePercentMatchingMotifs{c,4}=PercentMatchList{H};
                clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
            end
            
            if H>1;
                UniquePercentMatchingMotifs2=vertcat(UniquePercentMatchingMotifs2, UniquePercentMatchingMotifs);
                Uniquemotifmatchsizes2=vertcat(Uniquemotifmatchsizes2, Uniquemotifmatchsizes);
                ResortedCompiledStorage2=vertcat(ResortedCompiledStorage2, ResortedCompiledStorage);
            else
                UniquePercentMatchingMotifs2=UniquePercentMatchingMotifs;
                Uniquemotifmatchsizes2=Uniquemotifmatchsizes;
                ResortedCompiledStorage2=ResortedCompiledStorage;
            end
            clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts UniquePercentMatchingMotifs2 Uniquemotifmatchsizes2 ResortedCompiledStorage2 BlastDB ConcatenatedUniquePercentMatchingMotifs MaxNumMuts Reference Path_BlastDB
        catch
            H=H+2
            K=1;
            for y=1:size(TotalCompiledStorage,1)
                temp1=strsplit(TotalCompiledStorage{y,1},',');
                if str2double(temp1{3})==str2double(PercentMatchList{H})
                    temp2=strsplit(temp1{1,1},'_');
                    length(temp2{1,1});
                    qstart=str2num(temp1{7});
                    qend=str2num(temp1{8});
                    z=1;
                    for x=qstart:qend
                        temp3(z)=str2num(temp2{1,1}(x));
                        z=z+1;
                    end
                    PercentMatchingMotifs{K,1}=int2nt(temp3);
                    ResortedCompiledStorage(K,:)=TotalCompiledStorage(y,:);
                    motifmatchsizes{K,1}=temp1{4};
                    clearvars('temp3','temp2','temp1')
                    K=K+1;
                end
            end
            
            UniquePercentMatchingMotifs=PercentMatchingMotifs;
            Uniquemotifmatchsizes=motifmatchsizes;
            
            %Get Positions of Matching Motifs in Fasta File
            K=1;
            for Z=1:size(UniquePercentMatchingMotifs(:,:),1)
                M=0;
                V=1;
                while M==0
                    MotifMatchTemp=strfind(A(V).Sequence, PercentMatchingMotifs(Z));
                    if isempty(MotifMatchTemp)==0
                        MotifMatchCoordinates{K,1}=V;
                        MotifMatchCoordinates{K,2}=MotifMatchTemp;
                        MotifMatchCoordinates{K,3}=MotifMatchTemp+(length(UniquePercentMatchingMotifs{Z,1})-1);
                        K=K+1;
                        M=1;
                    end
                    V=V+1;
                end
                clearvars('M','V','z');
            end
            
            %Extract string of PolymorphismPos to compare next to PerfectMotifMatch
            for c=1:size(MotifMatchCoordinates(:,:),1)
                y=1;
                for x=MotifMatchCoordinates{c,2}:MotifMatchCoordinates{c,3};
                    MotifPolymorphismTemp(1,y)=polymorphismPos(MotifMatchCoordinates{c,1},x);
                    y=y+1;
                end
                MotifPolymorphismTemp=double(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=string(MotifPolymorphismTemp);
                MotifPolymorphismTemp2=strjoin(MotifPolymorphismTemp2);
                MotifPolymorphismTemp2= MotifPolymorphismTemp2{1}(find(~isspace(MotifPolymorphismTemp2)));
                UniquePercentMatchingMotifs{c,2}=MotifPolymorphismTemp2;
                UniquePercentMatchingMotifs{c,3}=sum(MotifPolymorphismTemp);
                UniquePercentMatchingMotifs{c,4}=PercentMatchList{H};
                clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
            end
            
            if H>1;
                UniquePercentMatchingMotifs2=vertcat(UniquePercentMatchingMotifs2, UniquePercentMatchingMotifs);
                Uniquemotifmatchsizes2=vertcat(Uniquemotifmatchsizes2, Uniquemotifmatchsizes);
                ResortedCompiledStorage2=vertcat(ResortedCompiledStorage2, ResortedCompiledStorage);
            else
                UniquePercentMatchingMotifs2=UniquePercentMatchingMotifs;
                Uniquemotifmatchsizes2=Uniquemotifmatchsizes;
                ResortedCompiledStorage2=ResortedCompiledStorage;
            end
            clearvars -except Path_BlastDB PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts UniquePercentMatchingMotifs2 Uniquemotifmatchsizes2 ResortedCompiledStorage2 BlastDB ConcatenatedUniquePercentMatchingMotifs MaxNumMuts Reference
        end
    end
end

UniquePercentMatchingMotifs2(:,5)=Uniquemotifmatchsizes2;

%Define Number of Mutations Explained by TRACE results
MotifsUsedByTRACE(:,1)=UniquePercentMatchingMotifs2(:,2);
MotifsUsedByTRACE(:,2)=UniquePercentMatchingMotifs2(:,1)
for x=1:size(MotifsUsedByTRACE,1)
    Temp={};
    for ggg=1:size(MotifsUsedByTRACE{x,1},2)
        Temp2=MotifsUsedByTRACE{x,1}(ggg);
        Temp=horzcat(Temp,Temp2);
    end
    Temp=transpose(str2num(str2mat(Temp)));
    for y=1:length(Temp)
        if Temp(y)==0
            Temp(y)=2;
        end
        if Temp(y)==1
            Temp(y)=3;
        end
    end
    TranslatedMotifsUsedByTRACE{x,1}=int2nt(Temp);
    TranslatedMotifsUsedByTRACE_Header{x,1}=sprintf('Motif_%d',x)
end

for x=1:size(MotifsUsedByTRACE,1)
    TranslatedMotifsUsedByTRACE0{x,1}=TranslatedMotifsUsedByTRACE{x,1};
    TranslatedMotifsUsedByTRACESeq{x,1}=MotifsUsedByTRACE{x,2};
end

fastawrite('MotifsUsedByTRACE_MutationLocation_Original.fas',TranslatedMotifsUsedByTRACE_Header,TranslatedMotifsUsedByTRACE0)
fastawrite('MotifsUsedByTRACE_Sequences_Original.fas',TranslatedMotifsUsedByTRACE_Header,TranslatedMotifsUsedByTRACESeq)
save('MotifsUsedByTRACE_Analysis_Original.mat')

A=fastaread('MotifsUsedByTRACE_Sequences_Original.fas');
ReferenceRegions={}
for x=1:size(A,1)
    x
    fastawrite('distIndex.fasta',A(x))
    [~,blastOutput]=system(sprintf('blastn.exe -db "%s%s" -query distIndex.fasta -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -word_size 11 -max_hsps 1 -max_target_seqs 1 -num_threads 2', Path_BlastDB, BlastDB));
    BlastTemp1=strsplit(blastOutput,',')
    delete distIndex.fasta
    pause(0.25)
    if strcmp(strtrim(BlastTemp1{13}),'plus');
        CoordinateStart=str2double(BlastTemp1{9});
        CoordinateEnd=str2num(BlastTemp1{10});
        for y=1:size(Reference,1)
            if strcmp(Reference(y).Header,BlastTemp1{2})==1
                k=1;
                for z=CoordinateStart:CoordinateEnd
                    ReferenceRegionsTemp(k)=Reference(y).Sequence(z);
                    k=k+1;
                end
                ReferenceRegionsTemp2{1}=(ReferenceRegionsTemp);
                ReferenceRegions{x,1}=ReferenceRegionsTemp2{1};
                clearvars ReferenceRegionsTemp ReferenceRegionsTemp2
            end
        end
    elseif strcmp(strtrim(BlastTemp1{13}),'minus');
        CoordinateStart=str2double(BlastTemp1{10});
        CoordinateEnd=str2num(BlastTemp1{9});
        for y=1:size(Reference,1)
            if strcmp(Reference(y).Header,BlastTemp1{2})==1
                k=1;
                for z=CoordinateStart:CoordinateEnd
                    ReferenceRegionsTemp(k)=Reference(y).Sequence(z);
                    k=k+1;
                end
                ReferenceRegionsTemp2{1}=seqrcomplement(ReferenceRegionsTemp);
                ReferenceRegions{x,1}=ReferenceRegionsTemp2{1};
                clearvars ReferenceRegionsTemp ReferenceRegionsTemp2
            end
        end
    end
    clearvars BlastTemp1
end
ReferenceRegions=upper(ReferenceRegions); %control for sequences that return lowercase
load('MotifsUsedByTRACE_Analysis_Original','MotifsUsedByTRACE')

MutationExplainedBlastMap={}
for x=1:size(MotifsUsedByTRACE)
    if length((ReferenceRegions{x}))-length(MotifsUsedByTRACE{x,2})~=0
        [Score,Alignment]=nwalign(ReferenceRegions{x},MotifsUsedByTRACE{x,2})
        ReferenceRegions{x}=Alignment(1,:);
        MotifsUsedByTRACE{x,2}=Alignment(3,:);
        if length(str2num(str2mat(MotifsUsedByTRACE{x,1})))-length(MotifsUsedByTRACE{x,2})~=0
            Intermediate={};
            for ggg=1:size(MotifsUsedByTRACE{x,1},2)
                Intermediate2=MotifsUsedByTRACE{x,1}(ggg);
                Intermediate=horzcat(Intermediate,Intermediate2);
            end
            Intermediate=transpose(str2num(str2mat(Intermediate)));
        end
    end
    testDataset=vertcat(ReferenceRegions{x},MotifsUsedByTRACE{x,2});
    polymorphismPos=testDataset~=repmat(testDataset(1,:), size(testDataset, 1), 1);
    polymorphismPos(testDataset=='-')=0;
    if strfind(ReferenceRegions{x},'-')~=0
        GapsToCorrect=strfind(ReferenceRegions{x},'-')
        for zz=1:length(GapsToCorrect)
            polymorphismPos(2,GapsToCorrect(zz))=0
        end
    end
    if strfind(MotifsUsedByTRACE{x,2},'-')~=0
        GapsToCorrect=strfind(MotifsUsedByTRACE{x,2},'-')
        for zz=1:length(GapsToCorrect)
            polymorphismPos(2,GapsToCorrect(zz))=0
            Intermediate=[Intermediate(:,1:GapsToCorrect(zz)-1) 0 Intermediate(:,GapsToCorrect(zz):end)];
            MotifsUsedByTRACE{x,1}=mat2str(Intermediate)
        end
    end
    
    %hotfix
    FixedMotif={};
    for ggg=1:size(MotifsUsedByTRACE{x,1},2)
        Intermediate2=MotifsUsedByTRACE{x,1}(ggg);
        FixedMotif=horzcat(FixedMotif,Intermediate2);
    end
    FixedMotif=transpose(str2num(str2mat(FixedMotif)));
    
    MutationExplainedBlastMap{x,1}=polymorphismPos(2,:)-FixedMotif;
    MutationExplainedPerHit{x,1}=sum(MutationExplainedBlastMap{x,1}(:) == -1);
    clearvars testDataset polymorphismPos Intermediate
end
UniquePercentMatchingMotifs2(:,6)=MutationExplainedPerHit

%Graph 2DHeatmap for Explained Muts
for x=1:size(UniquePercentMatchingMotifs2,1)
    ExpNumMuts(x)=UniquePercentMatchingMotifs2{x,6};
    PercentMatch(x)=str2num(UniquePercentMatchingMotifs2{x,4});
end

X=[ExpNumMuts;PercentMatch];
X=transpose(X);
hist3(X,'Ctrs',{0:1:MaxNumMuts 80:2:100},'CdataMode','auto')
xlabel('NumMuts')
ylabel('PercentMatch')
colorbar
view(2)

savefig(sprintf('%s_%s','HistogramHitsCorrected','Original'))


save('originaldata_by_matchlength.mat')
clear
%-----------------------------------------
%% Bin and Remove Original Data that are not above Background Set
load('originaldata_by_matchlength.mat')
originalmatchingCoordinates=cell2mat(UniquePercentMatchingMotifs2(:,6));
originalmatchingCoordinates(:,2)=cellfun(@str2double, UniquePercentMatchingMotifs2(:,4))
originalmatchingCoordinates(:,3)=cellfun(@str2double, UniquePercentMatchingMotifs2(:,5))
[o, oo, ooo]=unique(originalmatchingCoordinates, 'rows');

clearvars -except matchingCoordinates originalmatchingCoordinates o oo ooo
load('Backgrounddata_by_matchlength.mat')
modelmatchingCoordinates=cell2mat(UniquePercentMatchingMotifs2(:,6));
modelmatchingCoordinates(:,2)=cellfun(@str2double, UniquePercentMatchingMotifs2(:,4))
modelmatchingCoordinates(:,3)=cellfun(@str2double, UniquePercentMatchingMotifs2(:,5))
uniquemodelmatchingCoordinates=unique(modelmatchingCoordinates, 'rows');
[m, mm, mmm]=unique(modelmatchingCoordinates, 'rows');

originalfreq=accumarray(ooo,1)
originalfreq=originalfreq/sum(originalfreq);
originalmatchingfreq=originalfreq(ooo);
Backgroundfreq=accumarray(mmm, 1)
Backgroundfreq=Backgroundfreq/sum(Backgroundfreq);

for x=1:size(originalmatchingCoordinates,1);
    [check(x), pos(x)]=ismember(originalmatchingCoordinates(x,:), uniquemodelmatchingCoordinates, 'rows');
    if check(x)==1 && ((Backgroundfreq(pos(x))+2*std(Backgroundfreq))>originalmatchingfreq(x)) %%%% Adds 2 standard deviations for each Background data point, based on the entire Background data set
        cancel(x)=1;
    end
end

if exist('cancel')==1
    clearvars -except cancel
    
    load('originaldata_by_matchlength.mat')
    UniquePercentMatchingMotifs2(find(cancel==1), :)=[];
    ResortedCompiledStorage2(find(cancel==1),:)=[];
    
    FinalResultsMotif=UniquePercentMatchingMotifs2
    FinalResultsBlast=ResortedCompiledStorage2
    FinalResultsBlast(cell2mat(FinalResultsMotif(:,6))==0, :)=[];
    FinalResultsMotif(cell2mat(FinalResultsMotif(:,6))==0, :)=[];
    FinalResults=horzcat(FinalResultsBlast,FinalResultsMotif);
    save('FinalCleanedResults')
    
else
    load('originaldata_by_matchlength.mat')
    FinalResultsMotif=UniquePercentMatchingMotifs2
    FinalResultsBlast=ResortedCompiledStorage2
    FinalResultsBlast(cell2mat(FinalResultsMotif(:,6))==0, :)=[];
    FinalResultsMotif(cell2mat(FinalResultsMotif(:,6))==0, :)=[];
    FinalResults=horzcat(FinalResultsBlast,FinalResultsMotif);
    save('FinalCleanedResults')
end
end

