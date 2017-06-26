function parsave(fname, varargin)
    
    evalCommand = 'save(fname, ';
    for i = 2 : nargin
        evalCommand = [evalCommand '''' inputname(i) ''','];
        eval([inputname(i) '= varargin{i-1};']);
    end
    evalCommand = [evalCommand(1:end-1) ');'];
    eval(evalCommand);
    
end