function TRACE_IMGT_Formating(ExcelFile,ReferenceFasta,Clean_InputHandle)
% Script to take IMGT output and format for TRACE analysis. 
% Input Arguments:
% ExcelFile - tab delimited file from IMGT 2_IMGT-gapped-nt-sequences
% opened and saved as .xlsx in Microsoft Excel. 
% ReferenceFasta - IMGT reference IgHV file. Inlcuded as
% 'Human_IgHV_Reference_total.fasta'
% Clean_InputFileHandle - File extension for fasta input. Usually '*.fas'

%Generate Fastas from xls file
AutomatedFastaGen_IMGT_Gapped_NT(ExcelFile,ReferenceFasta);

%Clean fastas and append with identifier. Annotates Fasta files with
%'_forAnalysis.fas' if fasta file has >= 50 sequences. Removes gapped
%sequences and trims ends of sequence set. 
TRACE_FastaClean(Clean_InputHandle);

%Move files into respective folders
files = dir('*forAnalysis.fas')
for x=1:size(files,1)
    FolderName=files(x).name
    FolderName=erase(FolderName,'.fas_forAnalysis.fas')
    FolderName=FolderName(1:end-2)
    %clean Folder name
    for y=1:size(FolderName,1)
        FolderName=strtrim(FolderName)
        if contains(FolderName,'*')
            Index=strfind(FolderName,'*')
            for y=fliplr(Index:length(FolderName))
                FolderName(y)=[]
            end
        end
    end
    %Create Folder and place files inside
    mkdir(FolderName)
    copyfile('fastawriteBetter.m',FolderName)
    copyfile('TRACE_cycle.m',FolderName)
    copyfile('TRACE_noGaps_1000.m',FolderName)
    copyfile('TRACE_PostProcessing.m',FolderName)
    copyfile('ModelingMutationPos_CoordinateScramble.m',FolderName)
    copyfile('MutationCoverageAnalysis_PostTRACEv2.m',FolderName)
    copyfile('AutomatedAnalysisTRACE.m',FolderName)
    movefile(files(x).name,FolderName)
    %If more than one file per Folder, select largest file
    cd(FolderName)
    files2 = dir('*forAnalysis.fas')
    if size(files2,1)>1
        MaxFileSize=max(files2.bytes)
        for y=1:size(files2,1)
            if files2(y).bytes~=MaxFileSize
                delete(files2(y).name)
            end
        end
    end
    cd ..
end

%Remove unused fastas
delete *.fas
end