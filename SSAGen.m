function [time, Y] = SSAGen(SysInf, tfinal, recordStep, verbose_flag, split_flag)

numSpecies        = SysInf.numSpecies;
numReactions      = SysInf.numReactions;
speNames          = SysInf.speNames;
speValues         = SysInf.speValues;
cNames            = SysInf.cNames;
cValues           = SysInf.cValues;
speConstIndecies  = SysInf.speConstIndecies;
totalsIndecies    = SysInf.totalsIndecies;
VHolder           = SysInf.VHolder;

% initial values
X = speValues;

% extract V from VHolder and display
V = VHolder.V;

if verbose_flag
    disp( sprintf('Stoichiometric Matrix:\n') ); disp(V);
end

% matricies to hold propensity values, number of species present after each step, and the length of each step
A = zeros(numReactions,1);

% Actually do SSA ------------- %
[Y, X, time, run_time] = SingleTrajectory(V, X, speConstIndecies, numSpecies, speValues, tfinal, recordStep, verbose_flag);
% ----------------------------- %

if verbose_flag
	disp(sprintf('\nSpecies'' final amounts:\n'));
	Amount = X;
    dataTable = table(Amount,'RowNames',speNames);
    disp(dataTable);
end

if verbose_flag
    disp(' ');
end

end