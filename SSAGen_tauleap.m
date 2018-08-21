function [time, Y] = SSAGen(SysInf, tfinal, recordStep, verbose_flag, tau)

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

% attempt to pick tau if not specified - *extremely* crude estimate
if tau == 0
	% Single SSA trajectory to help determine good tau
	[Y, X, time, run_time] = SingleTrajectory(V, X, speConstIndecies, numSpecies, speValues, tfinal, recordStep, verbose_flag);
	tau = ( time(length(time)) / ( length(time)*recordStep ) ) * 3 ;
end

if verbose_flag
	disp( sprintf('Chosen value for tau:\n') ); disp(tau);
end

X = speValues;

tic
[Y, X, time, run_time] = SingleTrajectory_tauleap(V, X, speConstIndecies, numSpecies, speValues, tfinal, recordStep, verbose_flag, tau);
time_with_leap = toc;

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