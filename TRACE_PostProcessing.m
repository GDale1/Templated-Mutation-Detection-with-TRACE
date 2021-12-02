function TRACE_PostProcessing(species,Path_RequiredFiles)

if strcmp(species,'Human')
    Directory=sprintf('%s%s',Path_RequiredFiles,'Human Gene Coordinates Master List.txt')
    L=importdata(Directory);
elseif strcmp(species,'Mouse')
    Directory=sprintf('%s%s',Path_RequiredFiles,'Mouse Gene Coordinates Master List.txt')
    L=importdata(Directory);
end

for x=2:size(L, 1)
    N{x,:}=strsplit(L{x}, '"');
    N{x}(6)=[];
    N{x}(5)=[];
    N{x}(3)=[];
    findComma=find(N{x}{1}=='	');
    FullList{x,:}=strsplit(L{x}, '\t');
    N{x}{1}=N{x}{1}(findComma(2)+1:findComma(3)-1);
    O{x, 1}=N{x}(1);
    O{x, 2}=horzcat(transpose(strsplit(N{x}{2}, ',')),transpose(strsplit(N{x}{3}, ',')));
    O{x, 3}=FullList{x}{13};
    OO{x,1}=N{x}(1);
    OO{x,2}=horzcat(transpose(strsplit(N{x}{3}, ',')), transpose(strsplit(N{x}{2}, ',')));
end

%%%Create IndexLibrary to assign Gene and Exon/Intron Status
%%%Will Query based on Chr Identity and Corresponding 1/0 Indexing of O
disp('Creating ChrIndex')
ChrList=unique(vertcat(O{2:end,1}));
for y=1:length(ChrList)
ChrIndex{y,1}=ChrList{y};
ChrIndex{y,2}(1)=0;
for x=2:length(O)
if strcmp(string(ChrList{y}),string((O{x,1}{1})))==1
ChrIndex{y,2}(x)=1;
else
ChrIndex{y,2}(x)=0;
end
end
ChrIndex{y,2}=logical(ChrIndex{y,2});
end
disp('Fin ChrIndex')


O(1,:)=[];
OO(1,:)=[];
for x=1:size(OO,1)
    for y=1:size(OO{x,2},1)-1
        OOO{x,1}=OO{x,1};
        OOO{x,2}{y,1}=OO{x,2}{y,1};
        OOO{x,2}{y,2}=OO{x,2}{y+1,2};
    end
    OOO{x,3}=O{x,3};
end

disp('Begin Folder Cycling and Intron\Exon Assessment')

% Get a list of all files and folders in this folder.
files = dir
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% Print folder names to command window.
for z=3:size(subFolders, 1)
    disp(subFolders(z).name)
    clearvars storage
    cd(subFolders(z).name)
    matfile=dir('*.mat');
    if size(matfile, 1)>0
        load(matfile.name, 'storage');
        clearvars -except P O OOO L storage subFolders ChrIndex ChrList species
        if exist('storage', 'var')
            storage=horzcat(storage, repmat({[]}, size(storage,1), 1));
            for x=1:size(storage, 1)
                findComma=find(storage{x}==',');
                chromosome=storage{x}(findComma(1)+1:findComma(2)-1);
                coordinates(1)=str2double(storage{x}(findComma(8)+1:findComma(9)-1));
                coordinates(2)=str2double(storage{x}(findComma(9)+1:findComma(10)-1));
                storage{x,2}=[];
                for y=1:length(ChrList)
                    if strcmp(string(ChrList{y}),string(chromosome))==1
                       SpecificLines=find(ChrIndex{y,2});
                       for zz=1:length(SpecificLines)
                           for yy=1:size(O{SpecificLines(zz),2},1)
                                if (coordinates(1)>=str2double(O{SpecificLines(zz),2}{yy, 1}) && coordinates(1)<=str2double(O{SpecificLines(zz),2}{yy, 2}) || (coordinates(2)>=str2double(O{SpecificLines(zz),2}{yy, 1}) && coordinates(2)<=str2double(O{SpecificLines(zz),2}{yy, 2})))
                                storage{x,2}='exon';
                                storage{x,3}=O{SpecificLines(zz),3};
                                end
                           end
                           
                           if size(OOO{SpecificLines(zz),2},1)>1
                                for yyy=1:(size(OOO{SpecificLines(zz),2},1)-1)
                                     if strcmp(string(storage{x,2}),'exon')==1
                                            if (coordinates(1)>=str2double(OOO{SpecificLines(zz),2}{yyy, 1}) && coordinates(1)<=str2double(OOO{SpecificLines(zz),2}{yyy, 2}) || (coordinates(2)>=str2double(OOO{SpecificLines(zz),2}{yyy, 1}) && coordinates(2)<=str2double(OOO{SpecificLines(zz),2}{yyy, 2})))
                                            storage{x,2}='mixed';
                                            end
                                     elseif (coordinates(1)>=str2double(OOO{SpecificLines(zz),2}{yyy, 1}) && coordinates(1)<=str2double(OOO{SpecificLines(zz),2}{yyy, 2}) || (coordinates(2)>=str2double(OOO{SpecificLines(zz),2}{yyy, 1}) && coordinates(2)<=str2double(OOO{SpecificLines(zz),2}{yyy, 2})))
                                     storage{x,2}='intron';
                                     storage{x,3}=OOO{SpecificLines(zz),3};
                                     end
                                 end
                           end
                       end
                    end
                end
                if isempty(storage{x,2})==1
                   storage{x,2}='intergenic';
                   storage{x,3}=[];
                end
                
                temp=strsplit(storage{x,1},',');
                storage{x,4}=temp{13};
                storage{x,5}=temp{2};
                storage{x,6}=temp{4};
                clearvars temp
            end
            save('storage_PostProcessing.mat', 'storage');
        end
    end
    cd ..
end
disp('Finished Cycling')

disp('Begin Loading Processed Data')
CompiledStorage=[]
for z=3:size(subFolders, 1)
    cd(subFolders(z).name)
    matfile=dir('storage_PostProcessing.mat')
    
    if size(matfile, 1)>0
        load(matfile.name, 'storage');
        if min(size(storage))>0
        CompiledStorage=cat(1,CompiledStorage,storage)
        end
    end
   cd ..
end

%Percent Calcs
z=0
for x=1:size(CompiledStorage,1)
IntronCount=count(string(CompiledStorage{x,2}),'intron')
z=z+IntronCount
end
IntronCount=z
IntronPercent=IntronCount/size(CompiledStorage,1)

z=0
for x=1:size(CompiledStorage,1)
ExonCount=count(string(CompiledStorage{x,2}),'exon')
z=z+ExonCount
end
ExonCount=z
ExonPercent=ExonCount/size(CompiledStorage,1)

z=0
for x=1:size(CompiledStorage,1)
MixedCount=count(string(CompiledStorage{x,2}),'mixed')
z=z+MixedCount
end
MixedCount=z
MixedPercent=MixedCount/size(CompiledStorage,1)

z=0
for x=1:size(CompiledStorage,1)
IntergenicCount=count(string(CompiledStorage{x,2}),'intergenic')
z=z+IntergenicCount
end
IntergenicCount=z
IntergenicPercent=IntergenicCount/size(CompiledStorage,1)

z=0
for x=1:size(CompiledStorage,1)
PlusCount=count(string(CompiledStorage{x,4}),'plus')
z=z+PlusCount
end
PlusCount=z
PlusPercent=PlusCount/size(CompiledStorage,1)

z=0
for x=1:size(CompiledStorage,1)
MinusCount=count(string(CompiledStorage{x,4}),'minus')
z=z+MinusCount
end
MinusCount=z
MinusPercent=MinusCount/size(CompiledStorage,1)

%Chr Tally
v=[]
for x=1:size(CompiledStorage,1)
    t=string(CompiledStorage{x,5})
        if x==1
           v=t
        else
        v=vertcat(v,t)
        end
end
v=unique(v)
for k=1:length(v)
Chromosome_Dist{k,1}=v(k)
z=0
    for x=1:size(CompiledStorage,1)
    ChrCount=strcmp(string(CompiledStorage{x,5}),Chromosome_Dist{k,1})
    z=z+ChrCount
    end
Chromosome_Dist{k,2}=z
Chromosome_Dist{k,3}=z/(size(CompiledStorage,1))
end

%Lengths of Hits Tally
clearvars v ChrCount
v=[]
for x=1:size(CompiledStorage,1)
    t=string(CompiledStorage{x,6})
        if x==1
           v=t
        else
        v=vertcat(v,t)
        end
end
v=unique(v)
for k=1:length(v)
HitLength_Dist{k,1}=v(k)
z=0
    for x=1:size(CompiledStorage,1)
    ChrCount=count(string(CompiledStorage{x,6}),HitLength_Dist{k,1})
    z=z+ChrCount
    end
HitLength_Dist{k,2}=z
HitLength_Dist{k,3}=z/(size(CompiledStorage,1))
end


if exist('Chromosome_Dist')==0
Chromosome_Dist={}
end


if exist('HitLength_Dist')==0
HitLength_Dist={}
end

if exist('PlusPercent')==0
HitLength_Dist={}
end

if exist('MinusPercent')==0
HitLength_Dist={}
end


save('Compiled_PostProcessing.mat','CompiledStorage','IntronPercent','ExonPercent','MixedPercent','IntergenicPercent','Chromosome_Dist','HitLength_Dist','PlusPercent','MinusPercent')
 disp('Finished')
end
