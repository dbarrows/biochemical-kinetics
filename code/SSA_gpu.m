function Y = SSA_gpu(filename, SysInf, tfinal, verbose_flag)

numSpecies        = SysInf.numSpecies;
numReactions      = SysInf.numReactions;
speNames          = SysInf.speNames;
speValues         = SysInf.speValues;
cNames            = SysInf.cNames;
cValues           = SysInf.cValues;
speConstIndecies  = SysInf.speConstIndecies;
totalsIndecies    = SysInf.totalsIndecies;
VHolder           = SysInf.VHolder;
gpu = gpuDevice();

if verbose_flag
    disp( sprintf('GPU Device detected:\n') );
    disp( gpu );
end

V = VHolder.V;
vNumOfReactant = VHolder.vNumOfReactant;
vReactant = VHolder.vReactant;
vDimerMap = VHolder.vDimerMap;

fid = fopen('fireGpuTrajectories.m','w');
fprintf(fid,'function Y = fireGpuTrajectories(VHolder, verbose_flag)\n\n');

fprintf(fid, 'V = VHolder.V;\n\n');

for i = 1:numSpecies
    fprintf(fid, 'x%d = %d;\n', i, speValues(i));
end

fprintf(fid, '\ntfinal = %d;\n\n', tfinal);

fprintf(fid, '\tfunction [input');
for i = 1:numSpecies
    fprintf(fid, ', x%d', i);
end
fprintf(fid, ']');

fprintf(fid, ' = fire_single_gpu_trajectory(input');
for i = 1:numSpecies
    fprintf(fid, ', x%d', i);
end
fprintf(fid, ')\n\n');

format long

for i = 1:numReactions
    fprintf(fid, '\t\tc%d = %e;\n', i, cValues(i));
end

fprintf(fid, '\n\t\tt = 0;\n');

fprintf(fid, '\n\t\twhile t < tfinal\n\n');

for i = 1:numReactions

    % set initially to that reaction's parameter
    format long;
    fprintf(fid, '\t\t\ta%d = (%e)', i, cValues(i) );

    % multiply current value by each ractant's value if applicable, and account for dimerisation reactions
    for k = 1:(vNumOfReactant(i));
            
        curSpeciesIndex = vReactant(i,k);

        fprintf(fid,'*x%d', curSpeciesIndex);

        % determine is reactant is part of a dimerisation reaction using dimerisation map, then alter propensity accordingly
        dimer_number = vDimerMap(curSpeciesIndex, i);
        if dimer_number > 1
            count = 1;
            while (count < dimer_number)
                fprintf(fid,'*(x%d-%d)', curSpeciesIndex, count);
                count = count + 1;
            end

            fprintf(fid,'/%d', factorial(dimer_number) );
        end
    end

    fprintf(fid, ';\n', i, cValues(i) );
end

fprintf(fid, '\n\t\t\tasum =');
for i = 1:numReactions
    fprintf(fid, ' + a%d', i);
end
fprintf(fid, ';\n\n');

fprintf(fid, '\t\t\tif asum == 0\n\t\t\t\treturn\n\t\t\tend\n\n');

for i = 1:numReactions
    fprintf(fid, '\t\t\ta_tot_%d = (', i);
    for j = 1:i
        fprintf(fid,'+a%d', j);
    end
    fprintf(fid, ')/asum;\n');
end

fprintf(fid, '\n\t\t\tj = 1;\n\n');
fprintf(fid, '\t\t\trand_num = rand;\n\n');

fprintf(fid, '\t\t\tif a_tot_1 > rand_num\n');
fprintf(fid, '\t\t\t\tj = 1;\n');

for i = 2:numReactions
    fprintf(fid, '\t\t\telseif a_tot_%d > rand_num\n', i);
    fprintf(fid, '\t\t\t\tj = %d;\n', i);
end

fprintf(fid, '\t\t\tend\n\n');

fprintf(fid, '\t\t\ttau = log(1/rand)/asum;\n\n');

for i = 1:numSpecies
    if ~ismember(i,speConstIndecies)
        fprintf(fid, '\t\t\tx%d = x%d + V(%d,j);\n', i, i, i);
    end
end

fprintf(fid, '\n');
SBMLModel = TranslateSBML(filename);
count = 0;
for i = 1:length(SBMLModel.rule)
    curRule = SBMLModel.rule(i);

    for j = 1:length( speNames )
        if strcmp( speNames(j), curRule.variable )

            totalsIndecies(count+1) = j;
            count = count + 1;
            fprintf(fid, '\t\t\tx%d =', j);

            % token string array
            ruleToks = strsplit( curRule.formula ,'+');    

            for l = 1:length( ruleToks )
                for m = 1:length(speNames)
                    if strcmp( speNames(m), ruleToks(l) )
                        fprintf(fid, ' + x%d', m);
                    end
                end
            end

            fprintf(fid, ';\n', m);

            break;
        end
    end
end

fprintf(fid, '\n\t\t\tt = t + tau;\n\n');

fprintf(fid, '\t\tend\n\n');
fprintf(fid, '\tend');

fprintf(fid, '\n\nnum_trajectories = 10000;\n');

fprintf(fid, 'trial_nums = linspace(1,num_trajectories, num_trajectories)'';\n' );
fprintf(fid, 'inputs = gpuArray(trial_nums);\n' );

fprintf(fid, '\n[g_trial');
for i = 1:numSpecies
    fprintf(fid, ', g_x%d', i);
end
fprintf(fid, '] = arrayfun(@fire_single_gpu_trajectory, inputs');
for i = 1:numSpecies
    fprintf(fid, ', x%d', i);
end
fprintf(fid, ');\n\n');

fprintf(fid, 'trials = gather(g_trial);\n');

fprintf(fid, 'Y = [');
for i = 1:numSpecies
    fprintf(fid, ' gather(g_x%d)', i);
end
fprintf(fid, '];\n\n');

fprintf(fid,'\nend');

fclose(fid);

Y = fireGpuTrajectories(VHolder, verbose_flag);
wait(gpu);

Mean = zeros(numSpecies, 1);
Std_dev = zeros(numSpecies, 1);
for i = 1:numSpecies
    data = Y( : , i );
    Mean(i) = mean( data );
    Std_dev(i) = std( data );
end

if verbose_flag
    dataTableMean = table(Mean,'RowNames',speNames);
    disp(dataTableMean);
    dataTableStddev = table(Std_dev,'RowNames',speNames);
    disp(dataTableStddev);
end

end