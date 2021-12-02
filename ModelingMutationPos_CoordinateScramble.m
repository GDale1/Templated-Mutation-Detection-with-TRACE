function ModelingMutationPos_CoordinateScramble(InputFileHandle,inRange,NumMuts,NumReps)

names = dir(InputFileHandle);
names = {names.name};
Directory=erase(names,'.fas');

for cycle=1:NumReps
fprintf('%s %d\n', 'Cycle',cycle)
for cycle2=1:length(names)
disp (names{cycle2})
A=fastaread(names{cycle2});

% cd (Directory{cycle2});
% names2 = dir('*_TRACE*.mat');
% names2 = {names2.name};
% load(names2{1},'polymorphismPos');
% cd ..

testDataset=vertcat(A(:).Sequence);
polymorphismPos=testDataset~=repmat(testDataset(1,:), size(testDataset, 1), 1);
polymorphismPos(testDataset=='-')=0;

[rows,columns]=find(polymorphismPos);
UniqueRow=unique(rows);
CoordinateList=[rows,columns];
BackgroundCoordinateList=[0,0];
SlidingRecord=[0];
MaxRecord=[0];
NewSequenceList=char.empty(0,(size(polymorphismPos,2)));




    for x=1:length(UniqueRow)
    templist=[];
    Newtemplist=[];
  

        for xx=1:size(CoordinateList,1)
            if CoordinateList(xx,1)==UniqueRow(x)
            temp=CoordinateList(xx,2);
            templist=vertcat(templist,temp);
            end
        end
%templist contains all mutations coordiantes for a given sequence
        match=1;
        testDataset=vertcat(A(:).Sequence);
%polymorphismPos=testDataset~=repmat(testDataset(1,:), size(testDataset, 1), 1);
%polymorphismPos(testDataset=='-')=0;

        for y=1:size(testDataset, 2)-(inRange-1)
            if sum(polymorphismPos(UniqueRow(x), y:y+(inRange-1)))>=NumMuts
            temppositionStorage=find(polymorphismPos(UniqueRow(x), y:y+(inRange-1))==1);
            inrangepolymorphismPos{match, 1}=num2str(temppositionStorage+y-1);
            inrangepolymorphismPos_x(match,1)=UniqueRow(x);
            inrangepolymorphismPos_y(match,1)=y;
            inrangepolymorphismPos{match,2}=testDataset(UniqueRow(x), y:y+(inRange-1));
            match=match+1;
            end
        end

        if exist('temppositionStorage')==1
            [unique_inrangeMatches, index1, index2]=unique(inrangepolymorphismPos(inrangepolymorphismPos_x==UniqueRow(x),1));
        elseif exist('temppositionStorage')==0
            unique_inrangeMatches={};
        end

        if isempty(unique_inrangeMatches)==1
            MaxL=min(templist)-1;
            if MaxL<=0;
               MaxL==0;
            end
%max sliding for right movement
            MaxR=(size(polymorphismPos,2)-max(templist));
        elseif isempty(unique_inrangeMatches)==0
            for f=1:size(unique_inrangeMatches,1)
                for ff=1:size(((split(unique_inrangeMatches{f}))),1)
                    temp1=(split(unique_inrangeMatches{f}));
                    temp2(ff)=str2num(temp1{ff});
                end
                MaxInRangePositionList(f)=max(temp2);
                MinInRangePositionList(f)=min(temp2);
                clearvars temp1 temp2;
            end
            MaxR=(size(polymorphismPos,2))-max(MaxInRangePositionList);
            MaxL=min(MinInRangePositionList)-1;
        end

%define direction of sliding
        Dir=rand;
% Dir<=0.5 --> Negative Direction
% Dir>0.5 --> Positive Direction

%Uniformly slide mutation position
        if Dir<=0.5;
            SlidingDist=randi([0 MaxL]);
            Newtemplist=templist-SlidingDist;
            for kk=1:size(Newtemplist)
                if Newtemplist(kk)<=0
                    Newtemplist(kk)=0;
                end
            end

            for tt=1:(size(Newtemplist,1));
                if Newtemplist(tt)~=0
                   Newtemplist(tt,2)=UniqueRow(x);
                else
                   Newtemplist(tt,2)=0;
                end
            end

            SlidingRecord=vertcat(SlidingRecord,SlidingDist);
            MaxRecord=vertcat(MaxRecord,MaxL);

            for t=1:size(Newtemplist,1)
                ReverseNewtemplist(t,1)=Newtemplist(t,2);
                ReverseNewtemplist(t,2)=Newtemplist(t,1);
            end

            BackgroundCoordinateList=vertcat(BackgroundCoordinateList,ReverseNewtemplist);

        elseif Dir>0.5
            SlidingDist=randi([0 MaxR]);
            Newtemplist=templist+SlidingDist;
            for kk=1:size(Newtemplist)
                if Newtemplist(kk)>=size(polymorphismPos,2)
                   Newtemplist(kk)=0;
                end
            end

            for tt=1:(size(Newtemplist,1));
                if Newtemplist(tt)~=0
                   Newtemplist(tt,2)=UniqueRow(x);
                else
                   Newtemplist(tt,2)=0;
                end
            end

            SlidingRecord=vertcat(SlidingRecord,SlidingDist);
            MaxRecord=vertcat(MaxRecord,MaxR);

            for t=1:size(Newtemplist,1)
                ReverseNewtemplist(t,1)=Newtemplist(t,2);
                ReverseNewtemplist(t,2)=Newtemplist(t,1);
            end

            BackgroundCoordinateList=vertcat(BackgroundCoordinateList,ReverseNewtemplist);

        end

        Newtemplist2=Newtemplist(:,1);
        Newtemplist2=Newtemplist2(Newtemplist2~=0);
        %Pull reference fasta sequence and change mutations at each site in
        %Newtemplist
        clearvars m L
        m=char.empty;
        L=char.empty;
        %%%%%%%%%%%%%%%% Modded
        if Dir<=0.5;
            adjust=1*SlidingDist;
        else
            adjust=-1*SlidingDist;
        end
        
        for w=1:size(A(1).Sequence,2)-inRange+1
            if size(find(Newtemplist2>w & Newtemplist2<w+inRange-1),1)>=NumMuts
                checkRange(w)=1;
            end
        end
        if exist('checkRange', 'var');
                windows=find(checkRange==1);
                zz=1; 
                while zz<size(windows,2)
                windows((windows<=windows(zz)+inRange-1 & windows>windows(zz)))=[];
                zz=zz+1;
                end
%                 %clusterStart=find(checkRange(1:end-1)-checkRange(2:end)==-1);
%                 clusterStart=find(windows(1:end-1)-windows(2:end)==-1);
%                 clusterStart=clusterStart+1;
%                 for u=1:size(clusterStart,2)
%                     if size(clusterStart(clusterStart>=clusterStart(u) & clusterStart<clusterStart(u)+inRange-1))
%                         
%                     end
%                 end
clusterStart=windows;
                strstorage=A(1).Sequence;
                tempFind=find(Newtemplist(:,2)~=0);
                orig=A(Newtemplist(tempFind(1),2)).Sequence;
                NewSequence=A(1).Sequence;
                for w=1:size(clusterStart, 2)
                    if w~=size(clusterStart,2)
                        cluster{w}=windows((windows>=clusterStart(w) & windows<clusterStart(w+1)) & (windows+adjust+inRange-1 < length(A(1).Sequence)));
                    else
                        %cluster{w}=windows(windows>clusterStart(w) & windows<=windows(end));
                        cluster{w}=windows((windows>=clusterStart(w) & windows<=windows(end)) & (windows+adjust+inRange-1 < length(A(1).Sequence)));
                    end
                    if isempty(cluster{w})==0
                    whichWindow=ceil(rand(1)*size(cluster{w},2));
                    if cluster{w}(whichWindow)+adjust<=0 %enforce adjust results in real index in subsequent steps 
                       adjust=-1*cluster{w}(whichWindow)+1;
                       if (cluster{w}(whichWindow)+adjust)==0
                          adjust=adjust+1;
                       end
                    end
                    origmutations=find(orig(cluster{w}(whichWindow)+adjust:cluster{w}(whichWindow)+inRange-1+adjust)~=strstorage(cluster{w}(whichWindow)+adjust:cluster{w}(whichWindow)+inRange-1+adjust));
                    temporig=orig(cluster{w}(whichWindow)+adjust:cluster{w}(whichWindow)+inRange-1+adjust);
                    summuts=size(origmutations,2);
                    newmutations=setdiff(1:size(temporig,2), origmutations);
                    if summuts>size(newmutations, 2) %hotfix to catch error where randperm K>N
                       summuts=size(newmutations, 2)
                    end
                    newpos=newmutations(randperm(size(newmutations, 2), summuts));
                    muta=strstorage(cluster{w}(whichWindow):cluster{w}(whichWindow)+inRange-1);
                    for y=1:size(newpos, 2);
                        mutoptions(1:3,y)=setdiff([1 2 3 4], nt2int(muta(newpos(y))));
                        muta(1,newpos(y))= int2nt(mutoptions(ceil(rand(1)*3),y));
                    end
                    NewSequence(cluster{w}(whichWindow):cluster{w}(whichWindow)+inRange-1)=muta;
                    end
                end
        else
            strstorage=A(1).Sequence;
            k=repmat(['A','C','G','T'],size(Newtemplist2,1),1);
            for b=1:size(Newtemplist2,1)
                L(b,:)=k(b,k(b,:)~=strstorage(Newtemplist2(b)));
                m(1,b)=L(b,ceil(rand(1)*3));
            end
            
            NewSequence=A(1).Sequence;
            for z=1:length(m)
                NewSequence(Newtemplist2(z))=m(z);
            end
        end
        
        %%%%%%%%%%% End Modded
        %         strstorage=A(1).Sequence;
        %         k=repmat(['A','C','G','T'],size(Newtemplist2,1),1);
        %         for b=1:size(Newtemplist2,1)
        %             L(b,:)=k(b,k(b,:)~=strstorage(Newtemplist2(b)));
        %             m(1,b)=L(b,ceil(rand(1)*3));
        %         end
        
        %         NewSequence=A(1).Sequence;
        %         for z=1:length(m)
        %             NewSequence(Newtemplist2(z))=m(z);
        %         end
        
        NewSequenceList=vertcat(NewSequenceList,NewSequence);

        NewSeqs(x).Header=sprintf('%s%s', 'Background_',A(UniqueRow(x)).Header);
        NewSeqs(x).Sequence=NewSequence;
        clearvars unique_inrangeMatches Newtemplist Newtemplist2 MaxR MaxL ff SlidingDist Dir ReverseNewtemplist templist L k kk NewSequence inrangepolymorphismPos inrangepolymorphismPos_x inrangepolymorphismPos_y index1 index2 match temppositionStorage checkRange cluster mutoptions
        m=char.empty;
    end

    FinalSeqs=[A(1),NewSeqs];
    FinalSeq2=struct.empty;
    for t=1:size(FinalSeqs,2)
        if isempty(FinalSeqs(t).Header)==0
        FinalSeq2Add=FinalSeqs(t);
        FinalSeq2=[FinalSeq2, FinalSeq2Add];
        end
    end
    FinalSeqs=FinalSeq2;

    BackgroundCoordinateList(1,:)=[];
    SlidingRecord(1,:)=[];
    MaxRecord(1,:)=[];
    sortedCoordinateList=sortrows(CoordinateList,1);

fastawrite(sprintf('%s_%s',Directory{cycle2},'Background.fas'),FinalSeqs)
save(sprintf('%s_%s',Directory{cycle2},'Background.mat'))
clearvars -except names cycle Directory inRange NumMuts

end


mkdir(sprintf('%s_%d','Background',cycle))
movefile('*Background.fas',sprintf('%s_%d','Background',cycle))
movefile('*Background.mat',sprintf('%s_%d','Background',cycle))
cd (sprintf('%s_%d','Background',cycle))
mkdir('MAT Files')
movefile('*Background.mat','MAT Files')
cd ..
end
fprintf('%s\n', 'Finished')