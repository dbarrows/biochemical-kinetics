function [time, Y] = RREGen(SysInf, tfinal, verbose_flag, stiff_flag)

numSpecies        = SysInf.numSpecies;
numReactions      = SysInf.numReactions;
speNames          = SysInf.speNames;
speValues         = SysInf.speValues;
cNames            = SysInf.cNames;
cValues           = SysInf.cValues;
speConstIndecies  = SysInf.speConstIndecies;
totalsIndecies    = SysInf.totalsIndecies;
VHolder           = SysInf.VHolder;

if verbose_flag
    disp( sprintf('\nNumber of species types:\n') ); disp(numSpecies);
    disp( sprintf('\nNumber of reactions:\n') ); disp(numReactions);
end

if verbose_flag
    disp( sprintf('\nParameter names:\n') ); disp(cNames);
    disp( sprintf('\nParameter values:\n') ); disp(cValues);
end

if verbose_flag
    disp( sprintf('\nSpecies'' names:\n') ); disp(speNames);
    disp( sprintf('\nSpecies'' initial amounts:\n') ); disp(speValues);
end

% holder strings for RHS of RREs
for i = 1:numReactions
	A{i} = '';
end
A = A';

V = VHolder.V;
vNumOfReactant = VHolder.vNumOfReactant;
vReactant = VHolder.vReactant;
vDimerMap = VHolder.vDimerMap;

if verbose_flag
    disp( sprintf('Stoichiometric Matrix:\n') ); disp(V);
end

for i = 1:numReactions

    % set initially to that reaction's parameter
    format long;
    A{i} = strcat( A{i} , num2str( cValues(i) ) );
    
    % multiply current value by each reactant's value if applicable, and account for dimerisation reactions
    for k = 1:(vNumOfReactant(i));
        
        curSpeciesIndex = vReactant(i,k);

        A{i} = strcat( A{i} , sprintf('*X(%d)',curSpeciesIndex) );

        % determine if reactant is part of a dimerisation reaction using dimerisation map, then alter propensity accordingly
        dimer_number = vDimerMap(curSpeciesIndex, i);
        if dimer_number > 1
            for j = 1:(dimer_number-1)
                A{i} = strcat( A{i} , sprintf('*(X(%d)-%d)', curSpeciesIndex, j ) );
            end

            A{i} = strcat( A{i} , sprintf('/%d', factorial(dimer_number) ) );
        end
    end
end

fid = fopen('RRE_functions.m','w');
fprintf(fid,'function dXdt = RRE_functions(t,X)');

fprintf(fid,'\n\ndXdt = zeros(%d,1);\n\n', numSpecies);

for i = 1:numSpecies

	fprintf(fid,'dXdt(%d) = 0',i);

	if ~ismember(i,speConstIndecies)
		for j = 1:numReactions
			if V(i,j) ~= 0
				fprintf(fid, ' + %d*%s', V(i,j), A{j});
			end
		end
	end

	fprintf(fid,' ;\n');

end

fprintf(fid,'\nend');
fclose(fid);

%tspan = linspace(0,tfinal,tfinal);
tspan = [0 tfinal];

disp( sprintf('==========Starting Solver==========') )

if stiff_flag
    [time,y] = ode15s(@RRE_functions,tspan,speValues);
else
    [time,y] = ode45(@RRE_functions,tspan,speValues);
end

Y = y';

if length(totalsIndecies) ~= 0
	for i = 1:length(time)
		Y(:,i) = calculateSpecifiedTotals(Y(:,i));
	end
end

if verbose_flag
    disp(sprintf('\nSpecies'' final amounts:\n'));
    Amount = Y(:,length(Y));
    dataTable = table(Amount,'RowNames',speNames);
    disp(dataTable);
end

end