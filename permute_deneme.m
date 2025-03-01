clear

N = 5;
M = 6;

ns = (1:N)' + 5;
nps = 1:N;

ms = (1:M)' + 10;
mps = 1:M;

arg1 = (ns - nps);
arg1 = arg1(:);

arg2 = (ms - mps);
arg2 = arg2(:);
arg2 = arg2.';

C = arg1*arg2;

C = reshape(C, [N, N, M, M]); % ???

C2 = zeros([N, N, M, M]);

for n = ns'
    for np = nps
        for m = ms'
            for mp = mps
                temp = (n-np)*(m-mp);
                C2(n-5, np, m-10, mp) = temp;
            end
        end
    end
end

diff = C - C2; % ayni

A = rand(4,5,6,7);

b = (1:7)';
b = reshape(b, [1,1,1, 7])

size(A)
size(b)

A.*b;

s = sum(A,3);
size(s)