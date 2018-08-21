function SysInf = SSA_setup(filename, verbose_flag)

if verbose_flag
    disp(' ');
end

arg_type = class(filename);

switch arg_type
    case 'char'
        SBMLModel = TranslateSBML(filename);
    case 'struct'
        SBMLModel = filename;
end


% determine number of reactions and species present in the model
numReactions = length(SBMLModel.reaction);
numSpecies = length(SBMLModel.species);

if verbose_flag
    disp( sprintf('\nNumber of species types:\n') ); disp(numSpecies);
    disp( sprintf('\nNumber of reactions:\n') ); disp(numReactions);
end

[cNames, cValues] = GetParameters(SBMLModel);

if verbose_flag
    [cNames_us, cValues_us] = GetAllParameters(SBMLModel);
    cNames_us = cNames_us';
    cValues_us = cValues_us';

	disp(sprintf('\nParameter Values:\n'))
    Value = cValues_us;
    dataTable = table(Value,'RowNames',cNames_us);
    disp(dataTable);
end

% get species names and values, set X to initial values
[speNames, speValues] = GetSpecies(SBMLModel);
speNames = speNames';
speValues = speValues';

if verbose_flag
	disp(sprintf('\nSpecies'' initial amounts:\n'))
    Amount = speValues;
    dataTable = table(Amount,'RowNames',speNames);
    disp(dataTable);
end

% matricies to hold propensity values, number of species present after each step, and the length of each step
A = zeros(numReactions,1);

VHolder = StoichiometricMatricesHolder(SBMLModel);

% will generate 'calculatePropensities.m' file
GeneratePropensityCalculatorFile(SBMLModel, VHolder);

% will generate 'calculateSpecifiedTotals.m' file
totalsIndecies = GenerateSpecifiedTotalsCalculatorFile(SBMLModel);

% determine if any species have a boundary condition and get their indecies
speConstIndecies = GetConstantSpeciesIndecies(SBMLModel);

SysInf = SystemInformationHolder;

SysInf.numSpecies        = numSpecies;
SysInf.numReactions      = numReactions;
SysInf.speNames          = speNames;
SysInf.speValues         = speValues;
SysInf.cNames            = cNames;
SysInf.cValues           = cValues;
SysInf.speConstIndecies  = speConstIndecies;
SysInf.totalsIndecies    = totalsIndecies;
SysInf.VHolder           = VHolder;

end