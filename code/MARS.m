function varargout = MARS(filename, varargin)

% Options:
%
% 'Hist'    - generate a histogram of the results of 10,000 trajecories
% 'Verbose' - enable diagnostic outputs
% 'Time'    - max time to run the simulation for
% 'Record'  - step size between recording simulation state information
% 'Tau'     - use tau-leaping
% 'CLE'     - use Chemical Langevin Equation
% 'RRE'     - use raction rate equations
% 'GPU'     - use the system's CUDA-supported GPU for multiple trajectory generation
% 'MLMC'    - generate histogram using Multi-level Monte-Carlo simulation (experimental)
% 'Error'   - error to use for MLMC method
% 'Steps'   - number of data points to generate over the integration interval for all non-GPU parallel methods
% 'Graph'   - plot a graph of the results

% default values for user-provided arguments
hist_flag       = 0;
verbose_flag    = 0;
tau_flag        = 0;
cle_flag        = 0;
rre_flag        = 0;
stiff_flag      = 0;
split_flag      = 0;
keep_flag       = 0;
record_flag     = 0;
gpu_flag        = 0;
graph_flag      = 0;
mlmc_flag       = 0;
m_flag          = 0;
err_flag        = 0;
steps_flag      = 0;
tfinal          = 50;
recordStep      = 20;
numSteps        = 100;
err             = 0;
speciesToGraph  = [];

i = 1;
while (1+i) <= nargin
	switch varargin{i}
		case 'Hist'
			hist_flag = 1;

		case 'Verbose'
			verbose_flag = 1;

        case 'Time'
            if (i+2) > nargin
                disp( sprintf('\nSimulation time argument missing.\n') );
                return;
            else
                i = i + 1;
                time_arg = varargin{i};
                if ~isnumeric(time_arg)
                    disp( sprintf('\nSimulation time argument must be a number.\n') );
                    return;
                else
                    tfinal = time_arg;
                end
            end

        case 'Record' % get record step argument, check for validity (integer)
            if (i+2) > nargin
                disp( sprintf('\nRecord step size argument missing.\n') );
                return;
            else
                i = i + 1;
                record_arg = varargin{i};
                if ~isnumeric(record_arg) || mod(record_arg,1) ~= 0
                    disp( sprintf('\nRecord step size argument must be an integer.\n') );
                    return;
                else
                    recordStep = record_arg;
                    record_flag = 1;
                end
            end

        case 'Tau' % use tau-leaping, get value to use for tau
            if (i+2) > nargin
                tau = 0;
            else
                next_arg = varargin{i+1};
                if isnumeric(next_arg)
                    tau = next_arg;
                    i = i + 1;
                else
                    tau = 0;
                end
            end
            tau_flag = 1;

        case 'CLE' % use Langevin leaping algorithm, get value to use for tau
            if (i+2) > nargin
                tau = 0;
            else
                next_arg = varargin{i+1};
                if isnumeric(next_arg)
                    tau = next_arg;
                    i = i + 1;
                else
                    tau = 0;
                end
            end
            cle_flag = 1;

        case 'MLMC' % use MLMC method, get value to use for M, defaults to 100 steps
            if (i+2) > nargin
                M = 0;
            else
                m_arg = varargin{i+1};
                if isnumeric(m_arg) && mod(m_arg,1) ~= 0
                    disp( sprintf('\nMLMC M-value argument must be an integer.\n') );
                    return;
                elseif isnumeric(m_arg) && mod(m_arg,1) == 0
                    M = m_arg;
                    i = i + 1;
                    m_flag = 1;
                else
                    M = 0;
                end
            end
            mlmc_flag = 1;

        case 'Steps' % get specific numer of steps to generate for parallel methods not on a GPU
            if (i+2) > nargin
                disp( sprintf('\nNumber of steps size argument missing.\n') );
                return;
            else
                step_arg = varargin{i+1};
                if isnumeric(step_arg) && mod(step_arg,1) ~= 0
                    disp( sprintf('\nNumber of steps argument must be an integer.\n') );
                    return;
                elseif isnumeric(step_arg) && mod(step_arg,1) == 0
                    numSteps = step_arg;
                    i = i + 1;
                else
                    disp( sprintf('\nNumber of steps size argument missing or invalid.\n') );
                    return;
                end
            end

        case 'Error' % get error to use for MLMC method, overrides default
            if (i+2) > nargin
                disp( sprintf('\nError argument missing.\n') );
                return;
            else
                i = i + 1;
                err_arg = varargin{i};
                if ~isnumeric(err_arg)
                    disp( sprintf('\nError argument missing.\n') );
                    return;
                elseif isnumeric(err_arg) && err_arg < 0
                    disp( sprintf('\nError must be a positive number.\n') );
                else
                    err = err_arg;
                end
            end

        case 'Graph' % whether or not to grapht the results, and if so for which species (default is all)
            if (i+2) > nargin
                speciesToGraph = [];
            else
                next_arg = varargin{i+1};
                if isvector(next_arg) && isnumeric(next_arg)
                    speciesToGraph = next_arg;
                    i = i + 1;
                else
                    speciesToGraph = [];
                end
            end
            graph_flag = 1;

        case 'RRE'
            rre_flag = 1;

        case 'Stiff'
            stiff_flag = 1;

        case 'Split'
            split_flag = 1;

        case 'Keep'
            keep_flag = 1;

        case 'GPU'
            gpu_flag = 1;

		otherwise
			disp( sprintf('\nInvalid option detected at argument %d.\n',i) );
			return;
	end
	i = i + 1;
end
if (tau_flag + cle_flag + rre_flag + mlmc_flag + gpu_flag) > 1
    disp( sprintf('\nInvalid options: multiple methods selected. You may only pick one of CLE, Tau-Leaping, RRE, MLMC, or GPU\n') );
    return
end

if tfinal == 0
    disp( sprintf('\nInvalid final time argument\n') );
    return
end

if recordStep == 0
    disp( sprintf('\nInvalid record step size argument\n') );
    return
end

if numSteps == 0
    disp( sprintf('\nInvalid number of steps argument\n') );
    return
end

if m_flag && M < 2
    disp( sprintf('\nMLMC M-value must be an integer greater than 1\n') );
    return
end

addpath(genpath('./toolbox/SBMLToolbox'));

platform_str = computer;

% get type of platform so proper libraries can be added to path, currently only PC, Mac supported
switch platform_str
    case 'MACI64'
        addpath(genpath('./toolbox/libSBML/mac'));
    case 'PCWIN'
        addpath(genpath('./toolbox/libSBML/win32'));
    case 'PCWIN64'
        addpath(genpath('./toolbox/libSBML/win64'));
    case 'GLNXA64'
        addpath(genpath('./toolbox/libSBML/linux'));
    otherwise
        disp('Platform not supported');
        return;
end

% create files to be filled, then close all
file_list = {'calculatePropensities.m';...
            'calculateSpecifiedTotals.m';...
            'RRE_functions.m';...
            'fireGpuTrajectories.m'};

num_files = length(file_list);
for i = 1:num_files
    cur_file = file_list{i};
    fid = fopen(cur_file,'w');
    fclose(fid);
end

% get system (model) inforamtion from SBML file
SysInf = SSA_setup(filename, verbose_flag);
numSpecies        = SysInf.numSpecies;
numReactions      = SysInf.numReactions;
speNames          = SysInf.speNames;
speValues         = SysInf.speValues;
cNames            = SysInf.cNames;
cValues           = SysInf.cValues;
speConstIndecies  = SysInf.speConstIndecies;
totalsIndecies    = SysInf.totalsIndecies;
VHolder           = SysInf.VHolder;

% check for bad species-to-graph entries
if graph_flag && ~mlmc_flag
    numSpeToGraph = length(speciesToGraph);
    if numSpeToGraph ~= 0
        for i = 1:numSpeToGraph

            curIndex = speciesToGraph(i);
            if curIndex > numSpecies || curIndex < 1
                disp( sprintf(['Error: invalid species index provided, exceeds number of species in system or is less than 1.'...
                              ' Graph will not be displayed.\n']) );
                graph_flag = 0;
            end

            for j = 1:(i-1)
                checkIndex = speciesToGraph(j);
                if checkIndex == curIndex
                    disp( sprintf('Error: invalid species index provided, duplicate index. Graph will not be displayed.\n') );
                    graph_flag = 0;
                end
            end

        end
    end
end

method_name = '';
rehash

if gpu_flag
    Y = SSA_gpu(filename, SysInf, tfinal, verbose_flag);
    method_name = 'SSA on GPU';

elseif mlmc_flag
    %open parallel pool based on installed toolbox version
    version_less_flag = verLessThan('distcomp', '6.3');
    if version_less_flag
        matlabpool open;
    else
        parpool;
    end

    [Mean, Step] = MLMCGen(SysInf, tfinal, numSteps, verbose_flag, split_flag, speciesToGraph, graph_flag, M, err);
    Y = Mean;
    varargout{3} = Step;
    time = linspace(0,tfinal,numSteps);
    varargout{2} = time;
    method_name = 'MLMC';

    % close parallel pool
    if version_less_flag
        matlabpool close;
    else
        delete(gcp);
    end

elseif hist_flag
    %open parallel pool based on installed toolbox version
    version_less_flag = verLessThan('distcomp', '6.3');
    if version_less_flag
        matlabpool open;
    else
        parpool;
    end

    % use indicated method to generate trajectories accross multiple CPUs or CPU cores
    if tau_flag
        [Y,Mean,Std] = SSAGen_parfor_tauleap(SysInf, tfinal, recordStep, verbose_flag, tau, speciesToGraph, numSteps, graph_flag, 0);
        method_name = 'SSA with tau-leaping on parallel CPUs';
    elseif cle_flag
        [Y,Mean,Std] = SSAGen_parfor_cle(SysInf, tfinal, recordStep, verbose_flag, tau, speciesToGraph, numSteps, graph_flag);
        method_name = 'CLE on parallel CPUs';
    elseif rre_flag
        disp( sprintf('\n''Hist'' is not a valid option to use with the Reaction Rate Equation method\n') );
    else
	    [Y,Mean,Std] = SSAGen_parfor(SysInf, tfinal, recordStep, verbose_flag, speciesToGraph, numSteps, graph_flag);
        method_name = 'SSA on parallel CPUs';
    end

    varargout{2} = Mean;
    varargout{3} = Std;

    % close parallel pool
    if version_less_flag
        matlabpool close;
    else
        delete(gcp);
    end

else
    % generate single trajectory based on indicated method
    if tau_flag
        [time, Y] = SSAGen_tauleap(SysInf, tfinal, recordStep, verbose_flag, tau);
        method_name = 'SSA with tau-leaping';
    elseif cle_flag
        [time, Y] = SSAGen_cle(SysInf, tfinal, recordStep, verbose_flag, tau);
        method_name = 'CLE';
    elseif rre_flag
        if record_flag
            disp(sprintf('\nWarning: ''Record'' argument will be ignored - not valid with Reaction Rate Equation method\n'));
        end
        [time, Y] = RREGen(SysInf, tfinal, verbose_flag, stiff_flag);
        method_name = 'RRE';
    else
        [time, Y] = SSAGen(SysInf, tfinal, recordStep, verbose_flag);
        method_name = 'SSA';
    end

    varargout{2} = time;
end

varargout{1} = Y;

arg_type = class(filename);
switch arg_type
    case 'struct'
        filename = filename.name;
end

% graph end of simulation histogram is using a parallel method and graph flag has been set
if graph_flag
    if gpu_flag || hist_flag
        GraphHist(Y, speciesToGraph, speNames, split_flag, filename, method_name);
    else
        GraphPlot(Y, time, speciesToGraph, speNames, split_flag, filename, method_name);
    end
end

% keep temporary files if indicated
if ~keep_flag
    for i = 1:num_files
        cur_file = file_list{i};
        delete(cur_file);
    end
end

end