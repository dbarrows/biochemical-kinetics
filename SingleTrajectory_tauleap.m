function [Y, X, time, run_time] = SingleTrajectory_tauleap(V, X, speConstIndecies, numSpecies, speValues, tfinal, recordStep, verbose_flag, tau)

% set max muber of data points for each chunk of the recorded values matrix
numMaxDataPoints = 10000;

Y = zeros(numSpecies, numMaxDataPoints);
time = zeros(1, numMaxDataPoints);

% assign initial values recorded values
Y( : , 1 ) = X;
time(1) = 0;

% initial values
t = 0;
count = 1;

speConstIndeciesLength = length(speConstIndecies);

if verbose_flag
    disp( sprintf('\n=================STARTING SSA WITH TAU_LEAPING==============\n') );
end

tic

while t < tfinal

    %---------------------------------------------------------------
    % calculate value of each propensity at that step
    A = calculatePropensities(X);
    %---------------------------------------------------------------

    % break out of simulation if all species are consumed
    if cumsum(A) == 0
        if verbose_flag
            disp(sprintf('\n=========REACTION HALTED - ALL SPECIES COMSUMED========\n'));
            disp(sprintf('Final time was %f', t) );
        end

        Y( : , floor(count/recordStep) + 1) = X;
        time( floor(count/recordStep) + 1) = t;

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
    
    if mod(count, recordStep) == 0
        Y( : , (count/recordStep) + 1) = X;
        time(count/recordStep + 1) = t;
    end

    count = count + 1;
end

run_time = toc;

Y = Y( : ,1:(floor(count/recordStep)) );
time = time( 1:(floor(count/recordStep)) );

if verbose_flag
    disp( sprintf('\n%d steps taken\n', count-1 ) );
    disp( sprintf('%d steps recorded\n', floor(count/recordStep)-1 ) );
end

end