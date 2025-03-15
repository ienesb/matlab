function b = custom_mod(a, m)
    b = mod(a, m);
    if b == 0
        b = m;
    end
end