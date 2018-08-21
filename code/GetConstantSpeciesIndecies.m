function speConstIndecies = GetConstantSpeciesIndecies(SBMLModel)

numSpecies = length(SBMLModel.species);

speConstIndecies = zeros(numSpecies,1);

count = 0;
for i = 1:numSpecies
    if SBMLModel.species(i).constant
    	speConstIndecies(count + 1) = i;
    	count = count + 1;
    end
end

speConstIndecies = speConstIndecies(1:count,1);

end