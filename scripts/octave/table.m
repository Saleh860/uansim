classdef table
	properties (Access=public)
		Properties 
	end
	properties (GetAccess=public)
		data
	end
		
	methods (Access=public)
		function [this]=table(varargin)
			%constructor
			
			function [vars,varCount,varNames]=parseNameValuePairs(varargin)
				namedParameterIndex=find((cellfun(@ischar,varargin) & cellfun(@isrow,varargin)),1);
				if 	~isempty(namedParameterIndex)
					vars=varargin(1:namedParameterIndex-1);
					varCount=length(vars);
					p=inputParser();
					p.FunctionName='table';
					p.CaseSensitive=true;
					addParameter(p,'VariableNames',{}, validateNames(varCount));
					addParameter(p,'RowNames',{}, validateNames(rowCount));
					addParameter(p,'DimensionNames',{}, validateNames(2));
					p.KeepUnmatched=false;
					parse(p,varargin{namedParameterIndex:end});
					varNames=p.Results.VariableNames;
				else
					vars=varargin;
					varCount=length(vars);
					varNames={};		
				end
						
				if isempty(varNames)
					varNames=arrayfun(@(c)['Var', num2str(c)],1:varCount,'UniformOutput',false);
				end
			end
			
			function [v,d]=formatColumn(varName,col)
				[n,m]=size(col);
				if isa(col,'char') || (m==1)
					v={varName}; 
					if iscell(col)
						d=col;
					else 
						d=arrayfun(@(i){col(i,:)},(1:n)');
					end 
				else
					v={};
					d={};
					for j=1:m
						v{j}=[varName,'(',num2str(j),')'];
						if iscell(col)
							d(:,j)=col(:,j);
						elseif isnumeric(col)
							d(:,j)=arrayfun(@(i){col(i,j)},(1:n)');
						else
							error('Unexpected column type.')
						end
					end
				end
			end

			function vf=validateNames(varCount) 
				vf=@(vn)(iscell(vn) && size(vn,2)==varCount & all(cellfun(@ischar,vn) & cellfun(@isrow,vn)) );
            end

            this.Properties = TableProperties();
			rowCount=size(varargin{1},1);
			[vars,varCount,varNames]=parseNameValuePairs(varargin{:});
						
			n=size(vars{1},1);
			i=1;
			[this.Properties.VariableNames(i),this.data]=formatColumn(varNames{i}, vars{i});
			for i=2:varCount
				if size(vars{i},1)~=n
					error(['Error at column #', num2str(i),': All table columns must have the same number of rows'])
				else
					[v,d]=formatColumn(varNames{i}, vars{i});
					this.data=[this.data,d];
					this.Properties.VariableNames=[this.Properties.VariableNames,v];
				end
			end 
			
		end 
		
		function strTable=char(this)
			[n,m]=size(this.data);
			space=repmat(' ',n+2,2);
			for i=1:m
				var=this.data(:,i);
				if iscell(var)
					var=strvcat([this.Properties.VariableNames{i};cellfun(@num2str,var,'UniformOutput',false)]);
				elseif isnumeric(var)
					var=strvcat(this.Properties.VariableNames{i},num2str(var,'UniformOutput',false));
				elseif ischar(var)
					var=strvcat(this.Properties.VariableNames{i},var);
				else 
					error(['table::char: Unexpected column format: ',class(var)]);
				end				
				dashes=repmat('-',1,size(var,2));
				strCol=[var(1,:); dashes; var(2:end,:)];

				if i==1
					strTable=strCol;
					else
					strTable=[strTable,space,strCol];
				end
			end
		end
		
		function disp(this)
			disp(char(this))
		end
		
		function data=table2cell(this)
			data=this.data;
		end
	end
end
		
