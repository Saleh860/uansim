function b=myContains(str,pat)	if isOctave		if iscell(str)			b=~cellfun(@isempty,strfind(str,pat));		else 			b=~isempty(strfind(str,pat));		end
	else 		b=contains(str,pat);	end 	end 