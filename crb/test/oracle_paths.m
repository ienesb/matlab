function oracle_paths(action, which)
% ORACLE_PATHS  Add/remove legacy oracle directories for tests only.
%   oracle_paths('add', 'old')        old/eva/utils (63 verified functions)
%   oracle_paths('add', 'impl_otfs')  impl/otfs_ber
%   oracle_paths('add', 'impl_liu')   impl/2605.14451v1
%   oracle_paths('remove', ...)
% Production code never depends on these.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % matlab/
    switch which
        case 'old',       d = fullfile(root, 'old', 'eva', 'utils');
        case 'impl_otfs', d = fullfile(root, 'impl', 'otfs_ber');
        case 'impl_liu',  d = fullfile(root, 'impl', '2605.14451v1');
        otherwise, error('oracle_paths:which', 'unknown oracle %s', which);
    end
    assert(exist(d, 'dir') == 7, 'oracle_paths:missing', 'missing %s', d);
    switch action
        case 'add',    addpath(d);
        case 'remove', rmpath(d);
    end
end
