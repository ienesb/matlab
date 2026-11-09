function r = tc(r, name, err, tol)
% TC  Append one check (name, error, tolerance) to a results struct array.
    if isempty(r), r = struct('name', {}, 'err', {}, 'tol', {}, 'pass', {}); end
    k = numel(r) + 1;
    r(k).name = name;  r(k).err = err;  r(k).tol = tol;
    r(k).pass = isfinite(err) && (err <= tol);
    fprintf('  %-64s %11.3e  (tol %8.1e)  %s\n', name, err, tol, ternary(r(k).pass, 'ok', 'FAIL'));
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
