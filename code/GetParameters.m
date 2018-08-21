function [cNames, cValues] = GetParameters(SBMLModel)

numReactions = length(SBMLModel.reaction);

% get all parameter names, values
[allCNames, allCValues] = GetAllParameters(SBMLModel);

% create matricies to hold parameter names, values
cValues = zeros(numReactions,1);
cNames = cell(numReactions,1);

% get parameters for each individual reaction as they may not be in order
for i = 1:numReactions
    
    % attempt to get parameters from local rection context
    [cNamesTemp,cValuesTemp] = GetParameterFromReaction(SBMLModel.reaction(i));

    % if parameter values are NOT embedded in each local reaction context (return values will be NULL)
    if ( length(cNamesTemp) == 0 ) && ( length(cValuesTemp) == 0 )

        % declare parameter found flag, get string from reaction formula field, tokenize it by operator
        done = 0;
        kinLawString = SBMLModel.reaction(i).kineticLaw.formula;
        kinLawStringToks = strsplit( kinLawString ,'*');

        % searches for each token in the list of all parameters in the model, assigns those values to the correct reaction if found
        for j = 1:length(kinLawStringToks)
            curTok = kinLawStringToks(j);
            for k = 1:length(allCNames);
                if strcmp( allCNames{k} , curTok )
                    cNames{i} = allCNames{k};
                    cValues(i) = allCValues(k);
                    done = 1;
                    break;
                end

                if done == 1
                    break;
                end
            end

            if done == 1
                break;
            end
        end
    % if parameters were embedded in the local reaction contexts, assign them
    else
        cNames(i) = cNamesTemp;
        cValues(i) = cValuesTemp;
    end
end

end