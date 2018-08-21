function [Y, Step] = MLMCGen(SysInf, tfinal, numSteps, verbose_flag, split_flag, speciesToGraph, graph_flag, M, err)

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
numSteps = numSteps - 1;

if err == 0
	err = 1/100;
end

% extract V from VHolder and display
V = VHolder.V;

if verbose_flag
    disp( sprintf('Stoichiometric Matrix:\n') ); disp(V);
end

% matricies to hold propensity values, number of species present after each step, and the length of each step
A = zeros(numReactions,1);

% parameters
N = max(X);

% default M
if M == 0
	M = 3;
end

alp = zeros(numSpecies,1);
for i = 1:numSpecies
	val = X(i);
	if val == 0
		alp(i) = 0;
	else
		alp(i) = log(val)/log(N);
	end
end

bet = zeros(numReactions,1);
for i = 1:numReactions
	val = cValues(i);
	if val ~= 0
		bet(i) = log(val)/log(N);
	end
end

V_pos = -V;
V_pos(V_pos < 0) = 0;

% get gamma value from largest of candidates
gam = -Inf;
for i = 1:numSpecies
	for k = 1:numReactions
		if V(i,k) ~= 0
			gam_can = bet(k) + dot(V_pos(:,k),alp) - alp(i);
			if gam_can > gam
				gam = gam_can;
			end
		end
	end
end

% get rho value from largest of candidates
rho = Inf;
for k = 1:numReactions
	for i = 1:numSpecies
		if V(i,k) ~= 0
			rho_can = alp(i);
			if rho_can < rho
				rho = rho_can;
			end
		end
	end
end

L = ceil(abs(log(err)));

% should have three levels
if L <= 2
	l_0 = 0;
else
	l_0 = L - 2;
end

num_levels = L - l_0 + 1;

% get required number of coarsest trajectories
n_0 = 4 * ceil( (N^-rho * N^-gam * err^-2) / 4 );

% get level step sizes and required number of trajectories
h_l = zeros(num_levels, 1);
n_l = zeros(num_levels, 1);
for i = 1:num_levels
	l = l_0 + i - 1;
	h_l(i) = tfinal/(M^l);
	n_l(i) = ceil( N^-rho * N^gam * (L - l_0) * h_l(i) * err^-2 );
end

% make each n_l divisible by 4
for i = 1:num_levels
	val = n_l(i);
	while mod(val,4) ~= 0
		val = val+1;
	end
	n_l(i) = val;
end

% print information if required
if verbose_flag
	fprintf('N:\t%d\n', N);
	fprintf('Gamma:\t%d\n', gam);
	fprintf('Rho:\t%d\n', rho);
	fprintf('M:\t%d\n', M);
	fprintf('Error:\t%d\n', err);
	fprintf('Granularities:\n\n');
		disp(h_l);
	fprintf('Number of trajectories at each level:\n\n');
		disp(n_0);
		disp(n_l(2:num_levels));
end

num_trajectories = sum(n_l);

interval = tfinal/(numSteps+1);

% level 0
[~, Mean_coarse, ~] = SSAGen_parfor_tauleap(SysInf, tfinal, 0, 0, h_l(1), speciesToGraph, numSteps+1, 0, n_0);

Y = Mean_coarse;
time = linspace(0,tfinal,numSteps+1);

% ----------------------------- % start MLMC

for i = 2:num_levels

	num_runs = n_l(i);

	Y_sub = zeros( num_runs , numSpecies, numSteps+1);

	parfor k = 1:num_runs

		% setup that level trajectory
		hl 		= h_l(i);
		hl_1 	= M*hl;
		l 			= l_0 + i - 1;
		zl 			= X;
		zl_1 		= X;
		Zl			= zeros(numSpecies, numSteps+1);
		Zl_1		= zeros(numSpecies, numSteps+1);
		Zl(:,1) 	= X;
		Zl_1(:,1) 	= X;
		t 			= 0;
		n 			= 2;

		while n <= (numSteps+1)

			lam_bot = calculatePropensities(zl_1)';
			A = zeros(numReactions, 3);

			for j = 1:M

				lam_top = calculatePropensities(zl)';

				% (a)
				A(:,1) = min(lam_top, lam_bot);
				A(:,2) = lam_top - A(:,1);
				A(:,3) = lam_bot - A(:,1);

				% (b)
				Lam = poissrnd(A*hl);

				% (c)
				delta_top = Lam(:,1) + Lam(:,2);
				delta_bot = Lam(:,1) + Lam(:,3);
				zl = zl + V*delta_top;
				zl_1 = zl_1 + V*delta_bot;
                
                zl(speConstIndecies) = speValues(speConstIndecies);
                zl_1(speConstIndecies) = speValues(speConstIndecies);
                zl = calculateSpecifiedTotals(zl);
                zl_1 = calculateSpecifiedTotals(zl_1);

			end

			t = t + hl_1;

			% record
			if t > ((n-1)*interval)
				Zl(:,n) = zl;
				Zl_1(:,n) = zl_1;
				n = n + 1;
			end
		end
		
		data = Zl - Zl_1;
		Y_sub(k,:,:) = data;

	end

	Mean_level = zeros(numSpecies, numSteps+1);

	for i = 1:numSpecies
		for j = 1:(numSteps+1)
			Mean(i,j) = mean( Y_sub(:,i,j) );
		end 
	end

	Y = Y + Mean_level;

end

Step = h_l(length(h_l));

end