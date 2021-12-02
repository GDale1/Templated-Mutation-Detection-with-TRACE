function MutationCoverageAnalysis_PostTRACEv2(Type,InputFileHandle,MaxNumMuts)

%Load all necessary MAT data
if strcmp(Type,'Background')==1
    files = dir('Background*')
    dirFlags = [files.isdir]
    subFolders = files(dirFlags)
    
    for z=1:size(subFolders, 1)
        disp(subFolders(z).name)
        cd(subFolders(z).name)
        if exist('Compiled_PostProcessing.mat', 'file')==2
            matfile=dir('Compiled_PostProcessing.mat');
            if exist('CompiledStorage')==0
                load(matfile.name, 'CompiledStorage')
                TotalCompiledStorage=CompiledStorage
            else
                load(matfile.name, 'CompiledStorage')
                TotalCompiledStorage=vertcat(TotalCompiledStorage,CompiledStorage)
            end
        end
        cd ..
    end
end
%-----------------------------------------------------------------------
if strcmp(Type,'Original')==1
    files = dir('Original*')
    dirFlags = [files.isdir]
    subFolders = files(dirFlags)
    
    for z=1:size(subFolders, 1)
        disp(subFolders(z).name)
        cd(subFolders(z).name)
        if exist('Compiled_PostProcessing.mat', 'file')==2
            matfile=dir('Compiled_PostProcessing.mat');
            if exist('CompiledStorage')==0
                load(matfile.name, 'CompiledStorage')
                TotalCompiledStorage=CompiledStorage
            else
                load(matfile.name, 'CompiledStorage')
                TotalCompiledStorage=vertcat(TotalCompiledStorage,CompiledStorage)
            end
        end
        cd ..
    end
end
%Need Compiled FASTA for this to work
%------------------------------------------------------------------
%Generate Compiled Fasta
if strcmp(Type,'Background')==1
    % Get a list of all files and folders in this folder.
    files = dir('Background*')
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir]
    % Extract only those that are directories.
    subFolders = files(dirFlags)
    % Print folder names to command window.
    for z=1:size(subFolders, 1)
        disp(subFolders(z).name)
        cd(subFolders(z).name)
        names = dir(InputFileHandle)
        names = {names.name}
        for x=1:max(size(names))
            if exist('A')==0
                A=fastaread(names{x})
            else
                B=fastaread(names{x})
                A=vertcat(A,B)
            end
        end
        cd ..
    end
end
%----------------------------------------------------------------------
if strcmp(Type,'Original')==1
    % Get a list of all files and folders in this folder.
    files = dir('Original*')
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir]
    % Extract only those that are directories.
    subFolders = files(dirFlags)
    % Print folder names to command window.
    for z=1:size(subFolders, 1)
        disp(subFolders(z).name)
        cd(subFolders(z).name)
        names = dir(InputFileHandle)
        names = {names.name}
        for x=1:max(size(names))
            if exist('A')==0
                A=fastaread(names{x})
            else
                B=fastaread(names{x})
                A=vertcat(A,B)
            end
        end
        cd ..
    end
end
%End Data loading

%Generate Polymorphism Pos for Compiled FASTA
testDataset=vertcat(A(:).Sequence);
polymorphismPos=testDataset~=repmat(testDataset(1,:), size(testDataset, 1), 1);
polymorphismPos(testDataset=='-')=0;
polymorphismPos=double(polymorphismPos);


%Generate List of each level %match for each hit from
for y=1:size(TotalCompiledStorage,1)
    temp1=strsplit(TotalCompiledStorage{y,1},',');
    PercentMatchList{1,y}=temp1{3};
end
if exist('PercentMatchList')
    PercentMatchList=unique(cellfun(@str2num,PercentMatchList));
    PercentMatchList=cellfun(@string,(num2cell(PercentMatchList)));
    for x=1:length(PercentMatchList)
        if length(PercentMatchList{x})==2
            PercentMatchList{x}=strcat(PercentMatchList{x},'.000');
        end
        if strcmpi(PercentMatchList{x},'100')==1
            PercentMatchList{x}=strcat(PercentMatchList{x},'.000');
        end
        if length(PercentMatchList{x})==4
            PercentMatchList{x}=strcat(PercentMatchList{x},'00');
        end
        if length(PercentMatchList{x})==5
            PercentMatchList{x}=strcat(PercentMatchList{x},'0');
        end
    end
    
    %Get each level of Matching Motifs in NT form from TRACE output (each cycle
    %is for each percent and is stored as a plane of UniquePerfectMatchingMotifs
    for H=1:size(PercentMatchList,2)
        try
            K=1;
            for y=1:size(TotalCompiledStorage,1)
                temp1=strsplit(TotalCompiledStorage{y,1},',');
                if strcmp(temp1{3},PercentMatchList{H})==1
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
                    clearvars('temp3','temp2','temp1')
                    K=K+1;
                end
            end
            
            if exist('PercentMatchingMotifs')==1
                UniquePercentMatchingMotifs=unique(PercentMatchingMotifs,'stable');
                
                
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
                
                    %Keep rows in which Motif match coordinates are <=
                    %length of polymorphismPos
                    Index=find(cell2mat(MotifMatchCoordinates(:,3))<=size(polymorphismPos,2));
                    MotifMatchCoordinates=MotifMatchCoordinates(Index,:);
                    UniquePercentMatchingMotifs=UniquePercentMatchingMotifs(Index,:);
                
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
                    %clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2','PercentMatchingMotifs')
                    clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
                    
                end
                
                save(sprintf('%s_%s_%s%s','PerPercentMatchData',PercentMatchList{H},Type,'.mat'),'UniquePercentMatchingMotifs','MotifMatchCoordinates')
            end
            clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts
        catch
            try
                H=H+1
                K=1;
                for y=1:size(TotalCompiledStorage,1)
                    temp1=strsplit(TotalCompiledStorage{y,1},',');
                    if strcmp(temp1{3},PercentMatchList{H})==1
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
                        clearvars('temp3','temp2','temp1')
                        K=K+1;
                    end
                end
                
                if exist('PercentMatchingMotifs')==1
                    UniquePercentMatchingMotifs=unique(PercentMatchingMotifs,'stable');
                    
                    
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
                    
                    %Keep rows in which Motif match coordinates are <=
                    %length of polymorphismPos
                    Index=find(cell2mat(MotifMatchCoordinates(:,3))<=size(polymorphismPos,2));
                    MotifMatchCoordinates=MotifMatchCoordinates(Index,:);
                    UniquePercentMatchingMotifs=UniquePercentMatchingMotifs(Index,:);
                    
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
                        %clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2','PercentMatchingMotifs')
                        clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
                        
                    end
                    
                    save(sprintf('%s_%s_%s%s','PerPercentMatchData',PercentMatchList{H},Type,'.mat'),'UniquePercentMatchingMotifs','MotifMatchCoordinates')
                end
                clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts
            catch
                try
                if H>size(PercentMatchList,2)
                   H=find(cellfun(@str2num,PercentMatchList)==str2num(temp1{3}))
                end
                K=1;
                for y=1:size(TotalCompiledStorage,1)
                    temp1=strsplit(TotalCompiledStorage{y,1},',');
                    if strcmp(temp1{3},PercentMatchList{H})==1
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
                        clearvars('temp3','temp2','temp1')
                        K=K+1;
                    end
                end
                
                if exist('PercentMatchingMotifs')==1
                    UniquePercentMatchingMotifs=unique(PercentMatchingMotifs,'stable');
                    
                    
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
                    
                    %Keep rows in which Motif match coordinates are <=
                    %length of polymorphismPos
                    Index=find(cell2mat(MotifMatchCoordinates(:,3))<=size(polymorphismPos,2));
                    MotifMatchCoordinates=MotifMatchCoordinates(Index,:);
                    UniquePercentMatchingMotifs=UniquePercentMatchingMotifs(Index,:);
                    
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
                        %clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2','PercentMatchingMotifs')
                        clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
                        
                    end
                    
                    save(sprintf('%s_%s_%s%s','PerPercentMatchData',PercentMatchList{H},Type,'.mat'),'UniquePercentMatchingMotifs','MotifMatchCoordinates')
                end
                clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts
                catch
                if H>size(PercentMatchList,2)
                   H=find(cellfun(@str2num,PercentMatchList)==str2num(temp1{3}))
                end
                K=1;
                for y=1:size(TotalCompiledStorage,1)
                    temp1=strsplit(TotalCompiledStorage{y,1},',');
                    if strcmp(temp1{3},PercentMatchList{H})==1
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
                        clearvars('temp3','temp2','temp1')
                        K=K+1;
                    end
                end
                
                if exist('PercentMatchingMotifs')==1
                    UniquePercentMatchingMotifs=unique(PercentMatchingMotifs,'stable');
                    
                    
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
                    
                    %Keep rows in which Motif match coordinates are <=
                    %length of polymorphismPos
                    Index=find(cell2mat(MotifMatchCoordinates(:,3))<=size(polymorphismPos,2));
                    MotifMatchCoordinates=MotifMatchCoordinates(Index,:);
                    UniquePercentMatchingMotifs=UniquePercentMatchingMotifs(Index,:);
                    
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
                        %clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2','PercentMatchingMotifs')
                        clearvars('MotifPolymorphismTemp','MotifPolymorphismTemp2')
                        
                    end
                    
                    save(sprintf('%s_%s_%s%s','PerPercentMatchData',PercentMatchList{H},Type,'.mat'),'UniquePercentMatchingMotifs','MotifMatchCoordinates')
                end
                clearvars -except PercentMatchList Type TotalCompiledStorage A polymorphismPos MaxNumMuts
                
                
                end
            end 
        end
    end
    
    %load and concatenate individual files
    matfile2=dir(sprintf('%s%s%s','*',Type,'.mat'));
    for x=1:size(matfile2,1)
        load(matfile2(x).name);
        if exist('ConcatenatedUniquePercentMatchingMotifs')==0
            ConcatenatedUniquePercentMatchingMotifs=UniquePercentMatchingMotifs;
            ConcatenatedMotifMatchCoordinates=MotifMatchCoordinates;
        else ConcatenatedUniquePercentMatchingMotifs=vertcat(ConcatenatedUniquePercentMatchingMotifs,UniquePercentMatchingMotifs);
            ConcatenatedMotifMatchCoordinates=vertcat(ConcatenatedMotifMatchCoordinates,MotifMatchCoordinates);
        end
    end
    
    %Graph 2DHeatmap
    for x=1:size(ConcatenatedUniquePercentMatchingMotifs,1)
        NumMuts(x)=ConcatenatedUniquePercentMatchingMotifs{x,3};
        PercentMatch(x)=str2num(ConcatenatedUniquePercentMatchingMotifs{x,4});
    end
    
    X=[NumMuts;PercentMatch];
    X=transpose(X);
    hist3(X,'Ctrs',{0:1:MaxNumMuts 80:2:100},'CdataMode','auto')
    xlabel('NumMuts')
    ylabel('PercentMatch')
    colorbar
    view(2)
    
    savefig(sprintf('%s_%s','HistogramHits',Type))
    
    delete(sprintf('%s%s%s','*',Type,'.mat'))
    save(sprintf('%s%s%s','ConcatenatedResults_',Type,'.mat'))
end
end