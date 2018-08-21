function totalsIndecies =  GenerateSpecifiedTotalsCalculatorFile(SBMLModel)

numReactions = length(SBMLModel.reaction);
numSpecies = length(SBMLModel.species);

[speNames, speValues] = GetSpecies(SBMLModel);

totalsIndecies = zeros(numSpecies,1);

fid = fopen('calculateSpecifiedTotals.m','w');
fprintf(fid,'function X = calculateSpecifiedTotals(X)\n\n');

count = 0;
for i = 1:length(SBMLModel.rule)
	curRule = SBMLModel.rule(i);

    for j = 1:length( speNames )
        if strcmp( speNames(j), curRule.variable )

            totalsIndecies(count+1) = j;
            count = count + 1;
            fprintf(fid, 'X(%d) = 0', j);

            % token string array
            ruleToks = strsplit( curRule.formula ,'+');    

            for l = 1:length( ruleToks )
                for m = 1:length(speNames)
                    if strcmp( speNames(m), ruleToks(l) )
                        fprintf(fid, ' + X(%d)', m);
                    end
                end
            end

            fprintf(fid, ';\n', m);

            break;
        end
    end
end

fprintf(fid,'\nend',numReactions);

fclose(fid);

totalsIndecies = totalsIndecies(1:count);

end