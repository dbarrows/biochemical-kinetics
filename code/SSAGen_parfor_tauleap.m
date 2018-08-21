function varargout = SSAGen_parfor_tauleap(SysInf, tfinal, recordStep, verbose_flag, tau, speciesToGraph, numSteps, graph_flag, num_traj)

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
    disp(' ');
end

% set max muber of data points for each chunk of the recorded values matrix
if num_traj == 0
    numMaxDataPoints = 10008;
else
    numMaxDataPoints = num_traj;
end

% set X to initial values
X = speValues;

% matrix to hold propensity values, number of constant species
speConstIndeciesLength = length(speConstIndecies);
A = zeros(numReactions,1);

% extract V from VHolder and display
V = VHolder.V;
if verbose_flag
    disp( sprintf('\nStoichiometric Matrix:\n') ); disp(V);
end

if tau == 0
    % Single SSA trajectory to help determine good tau
    [Y, X, time, run_time] = SingleTrajectory(V, X, speConstIndecies, numSpecies, speValues, tfinal, recordStep, verbose_flag);
    tau = ( time(length(time)) / ( length(time)*recordStep ) ) * 3 ;
end

num_cores = feature('numCores');

if verbose_flag
    disp( sprintf('\nDetected %d CPU cores\n', num_cores) );
end

%Y = zeros(numMaxDataPoints, numSpecies);
Y = zeros(numSpecies, numSteps, numMaxDataPoints);
%time = zeros(1, numMaxDataPoints);

% begin SSA with tau-leaping algorithm

if verbose_flag
    disp( sprintf('\n============== STARTING parallel SSA with tau-leaping =============='  ) );
    disp( sprintf(  'Remember - this part takes a while. Please be patient.\n') );
end

halt_flag = 0;
time = linspace(0,tfinal,numSteps);
interval  = tfinal/numSteps;

parfor l = 1:numMaxDataPoints

    % initial values
    t = 0;
    n = 2;
    X = speValues;
    A = zeros(numReactions,1);
    Y_sub = zeros(numSpecies, numSteps);
    Y_sub(:,1) = X;

    while n <= numSteps

        next_step = (n-1)*interval;

        %---------------------------------------------------------------
        % calculate value of each propensity at that step
        A = calculatePropensities(X);
        %---------------------------------------------------------------

        asum = sum(A);

        % break out of simulation if all species are consumed
        if asum == 0
            halt_flag = 1;
            break
        end

        % get sampling of poisson random variables
        pois_rand_vars = poissrnd(A*tau);

        % update values
        X = X + V * pois_rand_vars';

        %---------------------------------------------------------------
        X = calculateSpecifiedTotals(X);
        %---------------------------------------------------------------

        for i = 1:speConstIndeciesLength
            X(speConstIndecies(i)) = speValues(speConstIndecies(i));
        end

        t = t + tau;

        if t > next_step
            Y_sub(:,n) = X';
            n = n + 1;
        end
        
    end

    Y(:,:,l) = Y_sub;
    
end

if halt_flag && verbose_flag
    disp(sprintf('\n=========REACTION HALTED - ALL SPECIES COMSUMED========\n'));
end

if verbose_flag
    disp(' ');
end

Y_last = zeros(numMaxDataPoints, numSpecies);

% get means, standard deviations, last slice data
Mean    = zeros(numSpecies, numSteps);
Std     = zeros(numSpecies, numSteps);
for i = 1:numSpecies
    for j = 1:numSteps
        data = Y(i,j,:);
        if j == numSteps
            Y_last(:,i) = data;
        end
        Mean(i,j) = mean(data);
        Std(i,j) = std(data);
    end
end

specific_species_flag = 0;

if ~isempty(speciesToGraph)
    stop_point = length(speciesToGraph);
    specific_species_flag  = 1;
else
    stop_point = numSpecies;
end
    
warning('off','MATLAB:legend:IgnoringExtraEntries');

if graph_flag
    for i = 1:stop_point
        if specific_species_flag
            index = speciesToGraph(i);
        else
            index = i;
        end
        
        data = zeros(numSteps, numMaxDataPoints);
        data(1:numSteps,1:numMaxDataPoints) = Y(index,:,:);
        data = data';
        
        figure;
        
        % graph boxplot for each slice
        subplot(2,1,1);
        boxplot( data , time, 'plotstyle', 'compact');
        legend( findobj(gca,'Tag','Box'),speNames(index) );
        xlabel('Time','FontSize',12, 'FontName', 'Helvetica');
        ylabel('Number of Species','FontSize',12,'FontName', 'Helvetica');
        title('Box and whisker plot of species vs time at given intervals','FontSize',16,'FontName', 'Helvetica');
        
        % graph means +- standard deviations for each slice
        subplot(2,1,2);
        hold all
        plot( time, Mean(index,:), 'b-o');
        plot( time, Mean(index,:) - Std(index,:), 'c');
        plot( time, Mean(index,:) + Std(index,:), 'c');
        legend( findobj(gca,'Tag','Box'),speNames(index) );
        xlabel('Time','FontSize',12, 'FontName', 'Helvetica');
        ylabel('Number of Species','FontSize',12,'FontName', 'Helvetica');
        title('Box and whisker plot of species vs time at given intervals','FontSize',16,'FontName', 'Helvetica');
    end
end

Mean_last = Mean(:,numSteps);
Std_dev_last = Std(:,numSteps);

varargout{1} = Y_last;
varargout{2} = Mean;
varargout{3} = Std;

if verbose_flag
    dataTableMean = table(Mean_last,'RowNames',speNames);
    disp(dataTableMean);
    dataTableStddev = table(Std_dev_last,'RowNames',speNames);
    disp(dataTableStddev);
end

end