function [strTable,table]=otable(varargin)
	namedParameterIndex=findFirstNamedParameter(varargin);
	if 	~isempty(namedParameterIndex)
		cols=varargin(1:namedParameterIndex-1);
		vCount=length(cols);
		p=inputParser();
		p.FunctionName='otable';
		p.CaseSensitive=true;
		addParameter(p,'VariableNames',{}, validateVariableNames(vCount));
		p.KeepUnmatched=false;
		parse(p,varargin{namedParameterIndex:end});
		vNames=p.Results.VariableNames;
	else
		cols=varargin;
		vCount=length(cols);
		vNames={};		
	end
	
	if isempty(vNames)
		vNames=arrayfun(@(c)(['Col #', num2str(c)]),1:vCount,'UniformOutput',false);
	end
	
	i=1;
	table=formatColumn(vNames{i}, cols{i});
	n=size(cols{i},1);
	for i=2:vCount
		if size(cols{i},1)~=n
			error(['Error at column #', num2str(i),': All table columns must have the same number of rows'])
		else
			v=formatColumn(vNames{i}, cols{i});
			table=[table,v];
		end
	end 
	m=size(table,2);
	colWidths=max([cellfun(@(c)(size(c,2)), table(2,:))
					cellfun(@length,table(1,:))]);
%	strTable=arrayfun(@(i)([table{1,i};table{2,i}]),1:vCount);
	i=1;
	dashes=repmat('-',1,colWidths(i));
	strTable=[table{1,i};dashes;table{2,i}];
	space=repmat(' ',n+2,2);
	for i=2:m
		dashes=repmat('-',1,colWidths(i));
		strTable=[strTable, space, [table{1,i};dashes;table{2,i}]];
	end
end

function v=formatColumn(heading,col)
	m=size(col,2);
	if m==1 || isa(col,'char')
		v={heading; num2str(col)};
	else
		v={};
		for i=1:m
			v=[v,{[heading,'(',num2str(i),')']; num2str(col(:,i))}];
		end
	end
end

function vf=validateVariableNames(varCount) 
	vf=@(vn)(iscell(vn) && size(vn,2)==varCount & all(cellfun(@ischar,vn) & cellfun(@isrow,vn)) );
end

function i=findFirstNamedParameter(v)
	i=find((cellfun(@ischar,v) & cellfun(@isrow,v)),1);
end
