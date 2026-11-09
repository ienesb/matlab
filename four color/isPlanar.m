function tf = isPlanar(A)
%ISPLANAR  Test whether an undirected graph is planar.
%   TF = ISPLANAR(A) returns logical true if the graph described by the
%   adjacency matrix A can be drawn in the plane with no edges crossing,
%   and false otherwise.
%
%   Input
%     A : N-by-N adjacency matrix. Only the nonzero pattern is used, so A
%         may be logical, 0/1, or weighted. A is treated as UNDIRECTED
%         (it is symmetrised internally). Self-loops and parallel edges are
%         ignored because they never affect planarity.
%
%   Method
%     1. Reduce A to a simple undirected graph.
%     2. Split the graph into its biconnected components (blocks). A graph
%        is planar iff every block is planar.
%     3. Test each block with the Demoucron-Malgrange-Pertuiset (DMP)
%        path-addition algorithm, after a quick Euler-bound rejection
%        (a simple planar graph on V>=3 vertices has at most 3V-6 edges).
%
%   The algorithm is exact (no heuristics): it returns the correct answer
%   for every graph, running in roughly O(V^2) time per block.
%
%   Examples
%     isPlanar(ones(5)-eye(5))                 % K5  -> false
%     isPlanar([zeros(3) ones(3); ones(3) zeros(3)])  % K3,3 -> false
%     isPlanar(ones(4)-eye(4))                 % K4  -> true

    A = (A ~= 0);
    if size(A,1) ~= size(A,2)
        error('isPlanar:square', 'Adjacency matrix must be square.');
    end
    A = A | A.';                 % treat as undirected
    n = size(A,1);
    A(1:n+1:end) = false;        % drop self-loops

    if n <= 4
        tf = true;               % every graph on <= 4 vertices is planar
        return;
    end

    adj = cell(n,1);
    for u = 1:n
        adj{u} = find(A(u,:));
    end

    blocks = biconnectedComponents(adj, n);

    tf = true;
    for b = 1:numel(blocks)
        if ~blockIsPlanar(blocks{b})
            tf = false;
            return;
        end
    end
end

% ======================================================================
function blocks = biconnectedComponents(adj, n)
% Iterative Hopcroft-Tarjan. Returns a cell array; each cell is a K-by-2
% edge list of one biconnected component (block).
    disc = zeros(1,n); low = zeros(1,n); parent = zeros(1,n);
    visited = false(1,n);
    timer = 0;
    estk = zeros(0,2);           % edge stack
    blocks = {};

    for s = 1:n
        if visited(s); continue; end
        fu = zeros(1,n); fi = zeros(1,n); top = 0;          % DFS frame stacks
        top = top+1; fu(top) = s; fi(top) = 1;
        visited(s) = true; timer = timer+1;
        disc(s) = timer; low(s) = timer; parent(s) = 0;

        while top > 0
            u = fu(top);
            if fi(top) <= numel(adj{u})
                v = adj{u}(fi(top));
                fi(top) = fi(top) + 1;
                if ~visited(v)
                    parent(v) = u;
                    estk(end+1,:) = [u v];                  %#ok<AGROW>
                    visited(v) = true; timer = timer+1;
                    disc(v) = timer; low(v) = timer;
                    top = top+1; fu(top) = v; fi(top) = 1;
                elseif v ~= parent(u) && disc(v) < disc(u)
                    estk(end+1,:) = [u v];                  %#ok<AGROW>
                    if disc(v) < low(u); low(u) = disc(v); end
                end
            else
                top = top - 1;
                if top > 0
                    p = fu(top);
                    if low(u) < low(p); low(p) = low(u); end
                    if low(u) >= disc(p)                    % p is a cut vertex / root
                        comp = zeros(0,2);
                        while true
                            e = estk(end,:); estk(end,:) = [];
                            comp(end+1,:) = e;              %#ok<AGROW>
                            if (e(1)==p && e(2)==u) || (e(1)==u && e(2)==p)
                                break;
                            end
                        end
                        blocks{end+1} = comp;               %#ok<AGROW>
                    end
                end
            end
        end
    end
end

% ======================================================================
function tf = blockIsPlanar(edges)
    vs = unique(edges(:)).';
    nv = numel(vs);
    ne = size(edges,1);
    if nv <= 4
        tf = true; return;
    end
    if ne > 3*nv - 6            % Euler bound: necessary condition for planarity
        tf = false; return;
    end
    remap = zeros(1, max(vs));  % relabel vertices to 1..nv
    remap(vs) = 1:nv;
    E = [remap(edges(:,1)).' remap(edges(:,2)).'];
    tf = dmp(E, nv);
end

% ======================================================================
function tf = dmp(E, nv)
% Demoucron-Malgrange-Pertuiset path-addition test on a biconnected graph.
    G = false(nv,nv);
    for i = 1:size(E,1)
        G(E(i,1),E(i,2)) = true; G(E(i,2),E(i,1)) = true;
    end
    adj = cell(nv,1);
    for u = 1:nv; adj{u} = find(G(u,:)); end

    H = false(nv,nv);           % embedded edges
    inEmb = false(1,nv);        % embedded vertices

    cyc = findCycle(adj, nv);   % seed embedding with any cycle
    k = numel(cyc);
    for i = 1:k
        a = cyc(i); b = cyc(mod(i,k)+1);
        H(a,b) = true; H(b,a) = true;
    end
    inEmb(cyc) = true;
    faces = {cyc, fliplr(cyc)};                 % interior + exterior face

    while true
        [fragEdges, contacts] = fragments(G, H, inEmb, nv);
        if isempty(fragEdges)
            tf = true; return;                  % everything embedded
        end

        nf = numel(faces);
        admCount = zeros(1, numel(fragEdges));
        admFirst = zeros(1, numel(fragEdges));
        for fi = 1:numel(fragEdges)
            c = contacts{fi};
            cnt = 0; first = 0;
            for fj = 1:nf
                if all(ismember(c, faces{fj}))
                    cnt = cnt + 1;
                    if first == 0; first = fj; end
                end
            end
            if cnt == 0
                tf = false; return;             % a fragment fits in no face
            end
            admCount(fi) = cnt; admFirst(fi) = first;
        end

        idx = find(admCount == 1, 1);           % DMP rule: forced fragment first
        if isempty(idx); idx = 1; end
        f = admFirst(idx);

        path = fragPath(fragEdges{idx}, contacts{idx}, inEmb, nv);
        faces = embedPath(faces, f, path);
        for i = 1:numel(path)-1
            H(path(i),path(i+1)) = true; H(path(i+1),path(i)) = true;
        end
        inEmb(path) = true;
    end
end

% ======================================================================
function [fragEdges, contacts] = fragments(G, H, inEmb, nv)
% Compute the bridges (fragments) of the embedded subgraph H within G.
    rem = G & ~H;
    comp = zeros(1,nv); cc = 0;
    for s = find(~inEmb)
        if comp(s) ~= 0; continue; end
        cc = cc + 1; q = s; comp(s) = cc; qi = 1;
        while qi <= numel(q)
            x = q(qi); qi = qi + 1;
            for y = find(rem(x,:))
                if ~inEmb(y) && comp(y) == 0
                    comp(y) = cc; q(end+1) = y;            %#ok<AGROW>
                end
            end
        end
    end

    compEdges = cell(1,cc); compContacts = cell(1,cc);
    fragEdges = {}; contacts = {};
    [ui, vi] = find(triu(rem));
    for t = 1:numel(ui)
        u = ui(t); v = vi(t);
        eu = inEmb(u); ev = inEmb(v);
        if ~eu && ~ev                                       % internal edge
            c = comp(u);
            compEdges{c}(end+1,:) = [u v];                  %#ok<AGROW>
        elseif eu && ev                                     % chord = own bridge
            fragEdges{end+1} = [u v];                       %#ok<AGROW>
            contacts{end+1} = [u v];                        %#ok<AGROW>
        else                                                % attachment edge
            if eu; w = u; c = comp(v); else; w = v; c = comp(u); end
            compEdges{c}(end+1,:) = [u v];                  %#ok<AGROW>
            compContacts{c}(end+1) = w;                     %#ok<AGROW>
        end
    end
    for c = 1:cc
        if isempty(compEdges{c}); continue; end
        fragEdges{end+1} = compEdges{c};                    %#ok<AGROW>
        contacts{end+1}  = unique(compContacts{c});         %#ok<AGROW>
    end
end

% ======================================================================
function path = fragPath(fe, ct, inEmb, nv)
% Find a path through one fragment connecting two of its contact vertices.
    if size(fe,1) == 1 && inEmb(fe(1,1)) && inEmb(fe(1,2))
        path = fe(1,:); return;                 % chord: the edge itself
    end
    adjF = cell(nv,1);
    for i = 1:size(fe,1)
        u = fe(i,1); v = fe(i,2);
        adjF{u}(end+1) = v; adjF{v}(end+1) = u;             %#ok<AGROW>
    end
    a = ct(1);
    pred = zeros(1,nv); seen = false(1,nv);
    q = a; seen(a) = true; qi = 1; target = 0;
    while qi <= numel(q)
        x = q(qi); qi = qi + 1;
        if x ~= a && inEmb(x)
            target = x; break;                  % reached another contact
        end
        if x == a || ~inEmb(x)                  % expand only through interior
            for y = adjF{x}
                if ~seen(y)
                    seen(y) = true; pred(y) = x; q(end+1) = y;  %#ok<AGROW>
                end
            end
        end
    end
    path = target;
    while path(1) ~= a
        path = [pred(path(1)) path];            %#ok<AGROW>
    end
end

% ======================================================================
function faces = embedPath(faces, f, path)
% Split face f along path, replacing it with the two resulting faces.
    F = faces{f};
    a = path(1); b = path(end);
    pa = find(F == a, 1); pb = find(F == b, 1);
    arc1 = forwardArc(F, pa, pb);               % a ... b
    arc2 = forwardArc(F, pb, pa);               % b ... a
    W = path(2:end-1);                          % new interior vertices
    faces(f) = [];
    faces{end+1} = [arc1 fliplr(W)];
    faces{end+1} = [arc2 W];
end

function arc = forwardArc(F, pi, pj)
    k = numel(F);
    if pi <= pj
        arc = F(pi:pj);
    else
        arc = [F(pi:k) F(1:pj)];
    end
end

% ======================================================================
function cyc = findCycle(adj, nv)
% Return the vertices of any cycle (assumes the graph contains one).
    visited = false(1,nv); parent = zeros(1,nv); disc = zeros(1,nv);
    st = zeros(1,nv); si = zeros(1,nv); top = 1;
    st(1) = 1; si(1) = 1; visited(1) = true; t = 1; disc(1) = t;
    while top > 0
        u = st(top);
        if si(top) <= numel(adj{u})
            v = adj{u}(si(top)); si(top) = si(top) + 1;
            if ~visited(v)
                visited(v) = true; t = t + 1; disc(v) = t; parent(v) = u;
                top = top + 1; st(top) = v; si(top) = 1;
            elseif v ~= parent(u) && disc(v) < disc(u)
                cyc = u; x = u;
                while x ~= v
                    x = parent(x); cyc(end+1) = x;          %#ok<AGROW>
                end
                return;
            end
        else
            top = top - 1;
        end
    end
    cyc = [];
end