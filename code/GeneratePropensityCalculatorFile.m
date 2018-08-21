% Generates a file 'calculatePropensities.m' that will contain a function able to calculate the necessary propensities for that system.

function GeneratePropensityCalculatorFile(SBMLModel, VHolder)

% determine number of reactions and species present in the model
numReactions = length(SBMLModel.reaction);
numSpecies = length(SBMLModel.species);

[cNames, cValues] = GetParameters(SBMLModel);

% matricies to hold propensity values
A = zeros(numReactions,1);

V = VHolder.V;
vNumOfReactant = VHolder.vNumOfReactant;
vReactant = VHolder.vReactant;
vDimerMap = VHolder.vDimerMap;

fid = fopen('calculatePropensities.m','w');
fprintf(fid,'function A = calculatePropensities(X)\n\n');

for i = 1:numReactions

    % set initially to that reaction's parameter
    format long;
    fprintf(fid, 'A(%d) = (%d)', i, cValues(i) );

    % multiply current value by each ractant's value if applicable, and account for dimerisation reactions
    for k = 1:(vNumOfReactant(i));
            
        curSpeciesIndex = vReactant(i,k);

        fprintf(fid,'*X(%d)', curSpeciesIndex);

        % determine is reactant is part of a dimerisation reaction using dimerisation map, then alter propensity accordingly
        dimer_number = vDimerMap(curSpeciesIndex, i);
        if dimer_number > 1
            count = 1;
            while (count < dimer_number)
                fprintf(fid,'*(X(%d)-%d)', curSpeciesIndex, count);
                count = count + 1;
            end

            fprintf(fid,'/%d', factorial(dimer_number) );
        end
    end

    fprintf(fid, ';\n', i, cValues(i) );
end

fprintf(fid,'\nend',numReactions);

fclose(fid);

end