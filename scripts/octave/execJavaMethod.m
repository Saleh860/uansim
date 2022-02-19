function result=execJavaMethod(className,methodName,varargin)	if isOctave		result=javaMethod(methodName,className,varargin{:});	else		f=eval(['@',className,'.',methodName]);		result=f(varargin{:});
	endend
