classdef StoichiometricMatricesHolder

    properties
        V;
        vNumOfReactant;
        vReactant;
        vDimerMap;
    end

    methods

        function object = StoichiometricMatricesHolder(SBMLModel)

            numReactions = length(SBMLModel.reaction);
            numSpecies = length(SBMLModel.species);

            % get the stoichiometric matrix representing the species changes for each reaction from the model
            V = zeros(numSpecies, numReactions);

            % matricies to hold the number of reactants in each reaction, the reactants' indecies, and the dimer map
            vNumOfReactant = zeros(numReactions,1);
            vReactant = zeros(numReactions, numSpecies);
            vDimerMap = zeros(numSpecies, numReactions);

            maxCount = 0;

            for j = 1:numReactions

                count = 0;
                for i = 1:numSpecies

                    role = DetermineSpeciesRoleInReaction(SBMLModel.species(i), SBMLModel.reaction(j));
                    if length(role) > 1

                        V(i,j) = role(1) - role(2);

                        if role(2) > 0
                            vReactant( j , (count+1) ) = i;
                            count = count + 1;
                        end

                        if role(2) >= 2
                            vDimerMap(i,j) = role(2);
                        end
                    else
                        V(i,j) = 0;
                    end
                end

                vNumOfReactant(j) = count;

                if count > maxCount
                    maxCount = count;
                end

            end
            
            object.V = V;
            object.vNumOfReactant = vNumOfReactant;

            % truncate reactant index matrix to discard unnecessary elements
            object.vReactant = vReactant( : , 1:(maxCount) );

            % make dimer reactant matrix sparse to discard unnecessary elements
            object.vDimerMap = vDimerMap;

        end

    end

end