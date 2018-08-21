function GraphHist(Y, speciesToGraph, speNames, split_flag, filename, method_name)


num_bins = length(Y)^(1/3);

specific_species_flag = 0;
speIndLength = length(speciesToGraph);

if speIndLength == 0
	end_point = min ( size(Y) );
else
	end_point = speIndLength;
	specific_species_flag = 1; 
end		

if ~split_flag
	figure
	hold all
end

for i = 1:end_point

	if specific_species_flag
		current_species = speciesToGraph(i);
	else
		current_species = i;
	end
		

	min_val = min( Y(:,current_species) );
	max_val = max( Y(:,current_species) );

	% determine optimal number of bins (educated guess)
    space = max_val - min_val + 1;
    if space > num_bins
        bins = linspace(min_val, max_val, num_bins);
        h = hist( Y(:,current_species) , num_bins );
    else
        bins = linspace(min_val, max_val, space);
        h = hist( Y(:,current_species) , space );
    end

	if split_flag
		figure
    end
    
	plot( bins, h/length(Y) , 'o-' );

	if split_flag
		legend( speNames(current_species) );
		xlabel('Number of Species','FontSize',12, 'FontName', 'Helvetica');
		ylabel('Frequency','FontSize',12,'FontName', 'Helvetica');
		title('Histogram of number of species at end of simulation','FontSize',16,'FontName', 'Helvetica');
	end
end

if ~split_flag
	hold off

	if specific_species_flag
		legend( speNames(speciesToGraph) );
	else
		legend( speNames )
	end

	xlabel('Number of Species','FontSize',12, 'FontName', 'Helvetica');
	ylabel('Frequency','FontSize',12,'FontName', 'Helvetica');

	title_string = ['Histogram of number of species at end of simulation from model source ''' filename ''' using ' method_name];
	title(title_string, 'Fontsize', 16,'FontName', 'Helvetica');
end

end