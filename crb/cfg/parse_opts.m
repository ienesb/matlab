function o = parse_opts(o, args, caller)
% PARSE_OPTS  Name/value override of a defaults struct that REJECTS unknown names
% (a misspelt option must not fall through silently).
%   o = parse_opts(o, varargin_cell, 'caller_name')
    if nargin < 3, caller = 'parse_opts'; end
    assert(mod(numel(args), 2) == 0, '%s:opts', 'options must come in name/value pairs', caller);
    for i = 1:2:numel(args)
        name = args{i};
        assert(ischar(name) || isstring(name), '%s:opts', 'option names must be strings', caller);
        name = char(name);
        assert(isfield(o, name), '%s:opt', '%s: unknown option ''%s'' (valid: %s)', caller, caller, name, strjoin(fieldnames(o).', ', '));
        o.(name) = args{i+1};
    end
end
