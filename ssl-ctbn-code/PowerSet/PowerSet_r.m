function [ P ] = PowerSet_r( S ,N)
    N
    x = 1:N;
    P = cell(1,2^N);
    p_ix = 2;
    for nn = 1:N
        a = combnk(x,nn);
        for j=1:size(a,1)
            P{p_ix} = S(a(j,:));
            p_ix = p_ix + 1;
        end
    end
end

