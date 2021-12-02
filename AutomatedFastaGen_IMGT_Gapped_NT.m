function AutomatedFastaGen_IMGT_Gapped_NT(ExcelFile,ReferenceFasta)
%For excel imported data using matlab popup
%OriginalVGeneList=Input(:,4);
%OriginalVGeneList=table2array(OriginalVGeneList);
%OriginalVGeneList=string(OriginalVGeneList);
[num,txt,Input]=xlsread(ExcelFile);

IndexPT1={};
for x=2:size(Input,1)
if isnan(Input{x,7})==1
IndexPT1{x-1}='Nan';
else
IndexPT1{x-1}=Input{x,7};
end
end
IndexPT1=string(IndexPT1);

for x=1:size(IndexPT1,2)
if strcmp(IndexPT1(x),'Nan')==1
IndexPT2(x)=0;
else
IndexPT2(x)=1;
end
end

y=1;
for x=1:size(IndexPT1,2)
    if IndexPT2(x)==1
    OriginalVGeneList{y}=string(Input{(x+1),4});
    IMGTgappedntsequences{y}=string(Input{(x+1),7});
    Header{y}=string(Input{(x+1),2});
    y=y+1;
    elseif IndexPT2(x)==0
    end
end
OriginalVGeneList=transpose(string(OriginalVGeneList));
IMGTgappedntsequences=transpose(string(IMGTgappedntsequences));
Header=transpose(string(Header));

Ref=fastaread(ReferenceFasta)

for x=1:(size(IMGTgappedntsequences,1))
if max(contains(OriginalVGeneList(x),', or'))==1
MutSeqsIndex(x)=0;
elseif max(contains(OriginalVGeneList(x),','))==0
MutSeqsIndex(x)=1;
end
end
%MutSeqsIndex(1)=[];
%OriginalVGeneList(1)=[];
UniqueIgHVMut=unique(OriginalVGeneList(MutSeqsIndex==1));

% Sequences=IMGTgappedntsequences(:,7);
% Sequences=table2array(Sequences);
% Sequences=string(Sequences);
% Sequences(1)=[];
Sequences=IMGTgappedntsequences;

%RefPreSplit{:}=Ref(:).Header;
for x=1:size(Ref,1)
RefPreSplit{1,x}=Ref(x).Header;
end
for z=1:length(RefPreSplit)
RefVGeneListTemp=strsplit((RefPreSplit{1,z}),'|');
RefVGeneList(z)=RefVGeneListTemp(2);
end

for y=1:length(UniqueIgHVMut)
%EachIgHv is indexed right to left down columns
for x=1:length(RefVGeneList)
Z=strsplit(UniqueIgHVMut(y));
EachIgHVIndex(y,x)=strcmp(RefVGeneList{x},Z(2));
end
end

disp 'Beginning Fasta Write'
%Create structures for fastas
for C=1:length(UniqueIgHVMut)
Alignment.Header=sprintf('%s_%s',erase(UniqueIgHVMut(C),'*'),'_Reference');
temp=EachIgHVIndex(C,:);
tempcheck=double(temp);
    if max(tempcheck)==1
        Alignment.Sequence=upper((Ref(temp==1).Sequence));
        for CC=1:length(Sequences)
            T=strsplit(OriginalVGeneList(CC));
            TT=strsplit(UniqueIgHVMut(C));
%if strcmp(T(2),TT(2))==1;
            if strcmp(erase(T(2),'*'),erase(TT(2),'*'))==1
             AlignmentMutated(CC).Sequence=upper(Sequences(CC));
             AlignmentMutated(CC).Header=Header(CC);
            end
        end
        if isempty(AlignmentMutated)==0
        Finalalignment=[Alignment,AlignmentMutated];
        end
        for j=1:length(Finalalignment)
        FinalIndex(j)=isempty(Finalalignment(j).Sequence);
        end
        Finalalignment2=Finalalignment(FinalIndex~=1)
        fastawrite(sprintf('%s%s',erase(TT(2),'*'),'.fas'),Finalalignment2)
        clearvars AlignmentMutated Alignment Finalalignment temp T TT Finalalignment2 FinalIndex
    end
end
disp 'End'
end

