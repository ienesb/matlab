function root = setup_crb()
% SETUP_CRB  Add the matlab/crb subfolders to the path and check for name
% collisions with the legacy library old/eva/utils.
%
% Production code never depends on old/eva; tests add it temporarily through
% test/oracle_paths.m.  Every function name in crb/ is chosen to differ from
% every old/eva/utils name, and this function asserts that no other copy of a
% crb/ function shadows or is shadowed on the current path.
%
% Usage:  cd matlab/crb; setup_crb; run_all_tests

    root = fileparts(mfilename('fullpath'));
    sub  = {'cfg','io','fim','rx','mcrb','util','run','plot','test'};
    for i = 1:numel(sub)
        d = fullfile(root, sub{i});
        if exist(d, 'dir'), addpath(d); end
    end

    % Name-collision guard: each crb/ function must resolve to exactly one file.
    for i = 1:numel(sub)
        d = fullfile(root, sub{i});
        if ~exist(d, 'dir'), continue; end
        files = dir(fullfile(d, '*.m'));
        for k = 1:numel(files)
            name = files(k).name(1:end-2);
            hits = which(name, '-all');
            hits = hits(~startsWith(hits, 'built-in') & ~contains(hits, 'is a variable'));
            if numel(hits) > 1
                error('setup_crb:collision', ...
                    'Function %s resolves to %d files on the path:\n%s', ...
                    name, numel(hits), strjoin(hits, newline));
            end
        end
    end
end
