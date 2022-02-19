function newStr = extractBefore(str,pat)	if ischar(pat)		if ischar(str)			if isrow(str)				k=strfind(str,pat);				if isempty(k)					newStr=str;				else 					newStr=str(1:k(1)-1);				end
			else				error('extractBefore: Invalid character input format.');			end
		elseif iscell(str)			error('extractBefore: Not yet implemented for cell arrays')		end
	else 		error('extractBefore: Not yet implemented for integer position')	endend 