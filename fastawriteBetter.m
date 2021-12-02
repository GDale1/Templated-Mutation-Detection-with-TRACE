        function fastawriteBetter(fileName, Header, Sequences);
            %tic
            formatSpec='%s\n';
            %tempHeader=Header;
            %tempSequences=Sequences;
            fileID=fopen(fileName, 'a');
            for x=1:size(Header, 1);
                fprintf(fileID, formatSpec, sprintf('>%s', Header(x,:)));
                fprintf(fileID, formatSpec, sprintf('%s', Sequences(x,:)));
            end
            fclose(fileID);
            %toc
        end