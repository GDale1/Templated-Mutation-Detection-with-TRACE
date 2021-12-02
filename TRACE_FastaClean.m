function TRACE_FastaClean(Clean_InputHandle)

names=dir(Clean_InputHandle)

for Cycle=1:size(names)
    A=fastaread(names(Cycle).name)
    
    %only proceed if file has >50 sequences
    if size(A,1)>=50
        
        %First minimize alignment
        k=1
        for x=1:length(A(1).Sequence)
            if strcmp(A(1).Sequence(x),'.')==1
                MinAlignIndex(k)=x
                k=k+1
            end
        end
        MinAlignIndex=fliplr(MinAlignIndex)
        
        for y=1:size(A,1)
            for x=1:length(MinAlignIndex)
                if strcmp(A(y).Sequence(MinAlignIndex(x)),'.')==1
                    A(y).Sequence(MinAlignIndex(x))=[];
                end
            end
        end
        
        %-----------------------------
        
        %Now remove CDR3
        for y=1:size(A,1)
            if y==1
                CutoffIndex=length(A(y).Sequence)
            end
            ToCut=CutoffIndex:length(A(y).Sequence);
            ToCut=fliplr(ToCut);
            for x=1:length(ToCut)
                if ToCut(x)>CutoffIndex
                    A(y).Sequence(ToCut(x))=[];
                end
            end
        end
        %-----------------------------
        
        %Now trim front end to remove gaps
        
        %Calculate the number of sequences with '.' as first character
        % NumberGappedSeqs=0
        % for y=1:size(A,1)
        %     if strcmp(A(y).Sequence(1),'.')==1
        %        NumberGappedSeqs=NumberGappedSeqs+1
        %     end
        % end
        
        InitialGapEndIndicies=[]
        for y=1:size(A,1)
            if strcmp(A(y).Sequence(1),'.')==1
                z=1
                while strcmp(A(y).Sequence(z),'.')==1
                    z=z+1;
                end
                InitialGapEndIndicies(y)=z;
            end
        end
        
        ChosenCutOff=median(InitialGapEndIndicies);
        CutOffMatrix=fliplr(1:ChosenCutOff);
        
        for y=1:size(A,1)
            for x=2:size(CutOffMatrix,2)
                A(y).Sequence(CutOffMatrix(x))=[];
            end
        end
        %---------------------------------------------
        
        %Now remove any sequence with a gap leftover
        CycleThrough=fliplr(1:size(A,1))
        for y=1:size(CycleThrough,2)
            if contains(A(CycleThrough(y)).Sequence,'.')
                A(CycleThrough(y))=[]
            end
        end
        %--------------------------------------------
        %only fastawrite if final fasta contains >=50 sequences
        if size(A,1)>=50
        fastawrite(sprintf('%s_%s',names(Cycle).name,'forAnalysis.fas'),A)
        end
        clearvars -except names Cycle
    end
end
end