function conditions = cartesianProduct(varargin)
    evalCommand = 'conditions = combvec(';
    for i = 1 : nargin
        evalCommand = [evalCommand '1:' inputname(i) ','];
        eval([inputname(i) '= varargin{i};']);
    end
    evalCommand = [evalCommand(1:end-1) ')'';'];
    eval(evalCommand);
