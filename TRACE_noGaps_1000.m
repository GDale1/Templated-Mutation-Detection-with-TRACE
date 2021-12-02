function TRACE_noGaps_1000(inputFasta, MotifLength, NumMuts, BlastDB, AnalysisName,Path_BlastDB)
A=fastaread(inputFasta);
testDataset=vertcat(A(:).Sequence);
polymorphismPos=testDataset~=repmat(testDataset(1,:), size(testDataset, 1), 1);
polymorphismPos(testDataset=='-')=0;
clearvars A;
tic
warning('off','all');
warning
inRange=MotifLength;
match=1;
storageIndex=1;
DuplicateList={};
for x=1:size(testDataset,1)
    if size(dir('*.fasta'),1)<=1000 & contains(inputFasta,'Background') %Impose upper threshold for number of motifs to check under Background
    for y=1:size(testDataset, 2)-(inRange-1)
        if sum(polymorphismPos(x, y:y+(inRange-1)))>=NumMuts
            temppositionStorage=find(polymorphismPos(x, y:y+(inRange-1))==1);
            inrangepolymorphismPos{match, 1}=num2str(temppositionStorage+y-1);
            inrangepolymorphismPos_x(match,1)=x;
            inrangepolymorphismPos_y(match,1)=y;
            inrangepolymorphismPos{match,2}=testDataset(x, y:y+(inRange-1));
            match=match+1;
        end
    end
    if exist('inrangepolymorphismPos_x', 'var')
        if any(inrangepolymorphismPos_x==x)
            [unique_inrangeMatches, index1, index2]=unique(inrangepolymorphismPos(inrangepolymorphismPos_x==x,1));
            for z=1:size(unique_inrangeMatches,1)
                subset=inrangepolymorphismPos(inrangepolymorphismPos_x==x,2);
                subsetx=inrangepolymorphismPos_x(inrangepolymorphismPos_x==x);
                subsety=inrangepolymorphismPos_y(inrangepolymorphismPos_x==x);
                string_fix=unique_inrangeMatches{z};
                string_fix(string_fix==' ')='_';
                fastaname=sprintf('%d_%s.fasta', x, string_fix);
                subsetIndex=find(index2==z);

                for headerIndex=1:size(subset(index2==z),1)
                    NT2NumHeader= join(num2str(nt2int(subset{subsetIndex(headerIndex)})));
                    NT2NumHeader= NT2NumHeader(find(~isspace(NT2NumHeader)));
                    fastaHeader{headerIndex}=sprintf('%s_%d_%d_%d', NT2NumHeader, x, subsety(subsetIndex(headerIndex)), subsetIndex(headerIndex));
                end
%% delete duplicate fastas
                subsetSubset=subset(subsetIndex);
                DuplicateAdd=vertcat(subsetSubset(:));                
                if isempty(DuplicateList)==0
                    for yy=1:length(subsetSubset)
                     if isempty(intersect(subsetSubset{yy},DuplicateList))==0
                        subsetSubset{yy}={};
                    end
                    end
                end
                
                for fastawriteindex=1:size(subsetSubset, 1)
                    if isempty(subsetSubset{fastawriteindex})==0
                    if sum(subsetSubset{fastawriteindex}=='-')<5
                       fastawriteBetter(fastaname, fastaHeader{fastawriteindex}, subsetSubset{fastawriteindex});
                    end
                    end
                end

                
                DuplicateList=vertcat(DuplicateList,DuplicateAdd);
                DuplicateList=unique(DuplicateList);

                if exist(fastaname)==2 
                    [~,blastOutput]=system(sprintf('blastn.exe -db "%s%s" -query %s -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -word_size 11 -max_hsps 1 -max_target_seqs 1 -num_threads 6', 'C:\Users\GDale\Desktop\Gordon\TRACE_Databases\Files\',BlastDB, fastaname));

                    if size(blastOutput,1)==0
                        delete(fastaname);
                    else
                        tempBlastOutput=splitlines(blastOutput);
                        for m=1:size(tempBlastOutput, 1)-1
                            if tempBlastOutput{m}(1)=='C'
                                toDelete(m)=m;
                            end
                        end
                        
                        if exist('toDelete', 'var')
                            tempBlastOutput(toDelete)=[];
                        end
                        
                        if size(tempBlastOutput, 1)>1
                        for y=1:size(tempBlastOutput, 1)-1
                            tempBlastcolumns{y}=strsplit(tempBlastOutput{y}, ',');
                        end
                        
                        tempBlastmat=vertcat(tempBlastcolumns{:});
                        sortedtempmat=sortrows(tempBlastmat, 12, 'descend');
                        topBlast=sortedtempmat(1,:);
                        storage2{storageIndex,1}=tempBlastOutput{1}; %%%%%%NOT CURRENTLY SELECTING THE HIGHEST BITSCORE!!!!
                        storage{storageIndex,1}=strjoin(topBlast, ',');
                        storageIndex=storageIndex+1;
                        end
                    end
                end
                clearvars fastaHeader subsetSubset tempBlastcolumns tempBlastmat sortedtempmat topBlast toDelete m
            end
        end
    end
    x
    elseif contains(inputFasta,'Background')==0
      for y=1:size(testDataset, 2)-(inRange-1)
        if sum(polymorphismPos(x, y:y+(inRange-1)))>=NumMuts
            temppositionStorage=find(polymorphismPos(x, y:y+(inRange-1))==1);
            inrangepolymorphismPos{match, 1}=num2str(temppositionStorage+y-1);
            inrangepolymorphismPos_x(match,1)=x;
            inrangepolymorphismPos_y(match,1)=y;
            inrangepolymorphismPos{match,2}=testDataset(x, y:y+(inRange-1));
            match=match+1;
        end
    end
    if exist('inrangepolymorphismPos_x', 'var')
        if any(inrangepolymorphismPos_x==x)
            [unique_inrangeMatches, index1, index2]=unique(inrangepolymorphismPos(inrangepolymorphismPos_x==x,1));
            for z=1:size(unique_inrangeMatches,1)
                subset=inrangepolymorphismPos(inrangepolymorphismPos_x==x,2);
                subsetx=inrangepolymorphismPos_x(inrangepolymorphismPos_x==x);
                subsety=inrangepolymorphismPos_y(inrangepolymorphismPos_x==x);
                string_fix=unique_inrangeMatches{z};
                string_fix(string_fix==' ')='_';
                fastaname=sprintf('%d_%s.fasta', x, string_fix);
                subsetIndex=find(index2==z);

                for headerIndex=1:size(subset(index2==z),1)
                    NT2NumHeader= join(num2str(nt2int(subset{subsetIndex(headerIndex)})));
                    NT2NumHeader= NT2NumHeader(find(~isspace(NT2NumHeader)));
                    fastaHeader{headerIndex}=sprintf('%s_%d_%d_%d', NT2NumHeader, x, subsety(subsetIndex(headerIndex)), subsetIndex(headerIndex));
                end
%%delete duplicate fastas
                subsetSubset=subset(subsetIndex);
                DuplicateAdd=vertcat(subsetSubset(:));
                
                if isempty(DuplicateList)==0
                    for yy=1:length(subsetSubset)
                     if isempty(intersect(subsetSubset{yy},DuplicateList))==0
                        subsetSubset{yy}={};
                    end
                    end
                end
                
                for xx=1:size(subsetSubset)
                if isempty(subsetSubset{xx})==0
                for fastawriteindex=1:size(subsetSubset, 1)
                    if isempty(subsetSubset{fastawriteindex})==0
                    if sum(subsetSubset{fastawriteindex}=='-')<5
                       fastawriteBetter(fastaname, fastaHeader{fastawriteindex}, subsetSubset{fastawriteindex});
                    end
                    end
                end
                end
                end
                
                DuplicateList=vertcat(DuplicateList,DuplicateAdd);
                DuplicateList=unique(DuplicateList);
                if exist(fastaname)==2 
                    [~,blastOutput]=system(sprintf('blastn.exe -db "%s%s" -query %s -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -word_size 11 -max_hsps 1 -max_target_seqs 1 -num_threads 6', Path_BlastDB,BlastDB, fastaname));

                    if size(blastOutput,1)==0
                        delete(fastaname);
                    else
                        tempBlastOutput=splitlines(blastOutput);
                        for m=1:size(tempBlastOutput, 1)-1
                            if tempBlastOutput{m}(1)=='C'
                                toDelete(m)=m;
                            end
                        end
                        
                        if exist('toDelete', 'var')
                            tempBlastOutput(toDelete)=[];
                        end
                        
                        if size(tempBlastOutput, 1)>1
                        for y=1:size(tempBlastOutput, 1)-1
                            tempBlastcolumns{y}=strsplit(tempBlastOutput{y}, ',');
                        end
                        
                        tempBlastmat=vertcat(tempBlastcolumns{:});
                        sortedtempmat=sortrows(tempBlastmat, 12, 'descend');
                        topBlast=sortedtempmat(1,:);
                        storage2{storageIndex,1}=tempBlastOutput{1}; %%%%%%NOT CURRENTLY SELECTING THE HIGHEST BITSCORE!!!!
                        storage{storageIndex,1}=strjoin(topBlast, ',');
                        storageIndex=storageIndex+1;
                        end
                    end
                end
                clearvars fastaHeader subsetSubset tempBlastcolumns tempBlastmat sortedtempmat topBlast toDelete m
            end
        end
    end  
    end
    x
end
toc

if exist('storage', 'var')>0
    
    
    for x=1:size(storage, 1)
        findUnderscore=find(storage{x}=='_');
        seqList{x}=(storage{x}(1:findUnderscore(1)-1));
        gapsList(x)=sum((seqList{x})=='6');
    end
    [shortList, indexA, indexB]=unique(seqList);
    
    neverCompressedStorage=storage;
    storage=storage(indexA);
    
    
    for y=1:size(storage, 1)
        bitscorePos=find(storage{y}==',');
        samplebitScore(y, 1)=str2num(storage{y}(bitscorePos(end-1)+1:bitscorePos(end)-1));
    end
    tic
    
    %% Mutate across window (Monte Carlo #1)
    for distIndex=1:size(storage, 1)
        distIndex
        
        tempfind=find(storage{distIndex}=='_');
        location=str2double(storage{distIndex}(tempfind(2)+1:tempfind(3)-1));
        locationStorage(distIndex)=location;
        sequenceNumber=str2double(storage{distIndex}(tempfind(1)+1:tempfind(2)-1));
        tempgermline=testDataset(1,location:location+inRange-1);
        tempgermline=tempgermline(testDataset(sequenceNumber, location:location+inRange-1)~='-');
        germlineSeq{distIndex,1}=tempgermline;
        for repetitions=1:1000
            germReps{distIndex, repetitions}=germlineSeq{distIndex,1};
            numMuts=sum(polymorphismPos(sequenceNumber,location:location+inRange-1));
            
            k=repmat(['A','C','G','T'],size(germlineSeq{distIndex},2),1);
            for i=1:1:size(germlineSeq{distIndex},2)
                L(i,:)=k(i,k(i,:)~=germlineSeq{distIndex}(i));
                m(1,i)=L(i,ceil(rand(1)*3));
            end
            
            mutSet=ceil(randperm(size(germReps{distIndex, repetitions}, 2), numMuts));
            germReps{distIndex, repetitions}(mutSet)=m(mutSet);
            
            NT2NumHeader2=join(num2str(nt2int(germReps{distIndex, repetitions})));
            NT2NumHeader2= NT2NumHeader2(find(~isspace(NT2NumHeader2)));
            
            fastawriteBetter('distIndex.fasta', sprintf('%s_%d',NT2NumHeader2,repetitions), germReps{distIndex, repetitions});
            
        end
        [~,distOutput{distIndex}]=system(sprintf('blastn.exe -db "%s%s" -query distIndex.fasta -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -word_size 11 -max_hsps 1 -num_threads 6', Path_BlastDB,BlastDB));
        delete distIndex.fasta
        toc
    end
    toc
    for distIndex=1:size(storage, 1)
        distIndex; 
        tempDist=distOutput{distIndex};
        tempDist(tempDist==char(10))=',';
        tempDist2=split(tempDist, ',');
        tempDist3=tempDist2(12:13:end);
        tempDist33=tempDist2(1:13:end);
        tempDist33(end)=[];
        
        tempDist333=split(tempDist33, '_');
        tempDist333=str2double(tempDist333);
        
        for x=1:size(tempDist3, 1)
            tempDist5(x)=str2num(tempDist3{x});
        end
        
        for x=1:max(tempDist333(:, 2))
            if size(max(tempDist5(tempDist333(:,2)==x)),2)>0
                tempDist4(x)=max(tempDist5(tempDist333(:,2)==x)); 
            end
        end
        tempDist4(tempDist4==0)=[];
        meanbitScore=mean(tempDist4);
        stdbitScore=round(std(tempDist4), 2);
        
        sampleZScore(distIndex)=(samplebitScore(distIndex, 1)-meanbitScore)/stdbitScore;
        clearvars tempDist5 tempDist4
    end
    sampleZScore_decompressed=sampleZScore(indexB);
    stoufferZ=sum(sampleZScore)/(size(sampleZScore,2)^0.5)
    stoufferZ_decompressed=sum(sampleZScore_decompressed)/(size(sampleZScore_decompressed,2)^0.5)
    
    clearvars -except storage germlineSeq sampleZScore stoufferZ polymorphismPos testDataset inRange BlastDB AnalysisName inputFasta Path_BlastDB
    originalStorage=storage;
    originalgermlineSeq=germlineSeq;
    originalsampleZScore=sampleZScore;
    originalstoufferZ=stoufferZ
    
    clearvars storage germlineSeq sampleZScore stoufferZ
    storage=originalStorage(originalsampleZScore>=1.645);
    germlineSeq=originalgermlineSeq(originalsampleZScore>=1.645);
    
    
    if isempty(germlineSeq)==0
        
        
        tic
        %% Mutate-In-Place (Monte Carlo #2)
        for distIndex=1:size(storage, 1)
            distIndex
            
            tempfind=find(storage{distIndex}=='_');
            location=str2double(storage{distIndex}(tempfind(2)+1:tempfind(3)-1));
            locationStorage(distIndex)=location;
            sequenceNumber=str2double(storage{distIndex}(tempfind(1)+1:tempfind(2)-1));
            tempgermline=testDataset(1,location:location+inRange-1);
            tempgermline=tempgermline(testDataset(sequenceNumber, location:location+inRange-1)~='-');
            germlineSeq{distIndex,1}=tempgermline;
            for repetitions=1:5
                germReps{distIndex, repetitions}=germlineSeq{distIndex,1};

                numMuts=sum(polymorphismPos(sequenceNumber,location:location+inRange-1)); 
                
                k=repmat(['A','C','G','T'],size(germlineSeq{distIndex},2),1);
                for i=1:1:size(germlineSeq{distIndex},2)
                    L(i,:)=k(i,k(i,:)~=germlineSeq{distIndex}(i));
                    m(1,i)=L(i,ceil(rand(1)*3));
                end

                germPos=nt2int(germlineSeq{distIndex,1});
                germPos=num2str(germPos);
                germPos(germPos==' ')=[];
                spacePos=find(storage{distIndex}=='_');
                tempStorage=storage{distIndex}(1:spacePos(1)-1);
                tempStorage=tempStorage(testDataset(sequenceNumber, location:location+inRange-1)~='-');
                mutSet=find(germPos~=tempStorage);
                germReps{distIndex, repetitions}(mutSet)=m(mutSet);
                

                NT2NumHeader2=join(num2str(nt2int(germReps{distIndex, repetitions})));
                NT2NumHeader2= NT2NumHeader2(find(~isspace(NT2NumHeader2))); 
                
                fastawriteBetter('distIndex.fasta', sprintf('%s_%d',NT2NumHeader2,repetitions), germReps{distIndex, repetitions});

            end
            [~,distOutput{distIndex}]=system(sprintf('blastn.exe -db "%s%s" -query distIndex.fasta -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -word_size 11 -max_hsps 1 -num_threads 6', Path_BlastDB,BlastDB));
            delete distIndex.fasta
            toc
        end
        toc
        toc
        for y=1:size(storage, 1)
            bitscorePos=find(storage{y}==',');
            samplebitScore(y, 1)=str2num(storage{y}(bitscorePos(end-1)+1:bitscorePos(end)-1));
        end
        tic
        
        
        
        for distIndex=1:size(storage, 1)
            distIndex;
            if size(distOutput{distIndex},1)>0
            tempDist=distOutput{distIndex};
            tempDist(tempDist==char(10))=','; 
            tempDist2=split(tempDist, ',');
            tempDist3=tempDist2(12:13:end);
            tempDist33=tempDist2(1:13:end);
            tempDist33(end)=[];
            tempDist333=split(tempDist33, '_');
            tempDist333=str2double(tempDist333);
            
            for x=1:size(tempDist3, 1)
                tempDist5(x)=str2num(tempDist3{x}); 
            end
            
            if size(tempDist333,1)==2 && size(tempDist333,2)==1
                tempDist333(1,2)=tempDist333(2,1);
                tempDist333(2,:)=[];
            end
            
            for x=1:max(tempDist333(:, 2))
                if size(max(tempDist5(tempDist333(:,2)==x)),2)>0
                    tempDist4(x)=max(tempDist5(tempDist333(:,2)==x));
                end
            end
            tempDist4(tempDist4==0)=[];

            meanbitScore=mean(tempDist4);
            stdbitScore=round(std(tempDist4), 2);
            
            sampleZScore(distIndex)=(samplebitScore(distIndex, 1)-meanbitScore)/stdbitScore;
            clearvars tempDist4 tempDist5
            end
        end
        stoufferZ=sum(sampleZScore)/((size(sampleZScore,2))^0.5)
    else
    end
else
end

save(sprintf('%s_TRACE_%s.mat', AnalysisName, inputFasta));
runtime=toc;
save('runtime.mat', 'runtime');

end
        