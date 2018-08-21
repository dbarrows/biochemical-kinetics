function GraphResults(Y, time, speciesToGraph, speNames, split_flag, filename, method_name)

specific_species_flag = 0;
speIndLength = length(speciesToGraph);

if speIndLength == 0
	dims = size(Y);
	end_point = dims(1);
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
		
	if split_flag
		figure
	end
		
	plot( time , Y(current_species,:) );

	if split_flag
		legend( speNames(current_species) );
		xlabel('Time','FontSize',12, 'FontName', 'Helvetica');
		ylabel('Number of Species','FontSize',12,'FontName', 'Helvetica');
		title('Species vs time using SSA','FontSize',16,'FontName', 'Helvetica');
	end
end

if ~split_flag
	hold off
	
	if specific_species_flag
		legend( speNames(speciesToGraph) );
	else
		legend( speNames )
	end

	xlabel('Time','FontSize',12, 'FontName', 'Helvetica');
	ylabel('Number of Species','FontSize',12,'FontName', 'Helvetica');

	title_string = ['Species vs time from model source ''' filename ''' using ' method_name];
	title(title_string,'FontSize',16,'FontName', 'Helvetica');
end

end