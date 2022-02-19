function newStr=myJoin(str,delimiter)	if ~isOctave		newStr=cell2mat(join(str,delimiter));	elseif iscell(str) && ischar(delimiter)		newStr=[str{1}];		for i=2:length(str)			newStr=[newStr,delimiter,str{i}];		end
	end
end 