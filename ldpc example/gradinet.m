xs = 5;
alpha = 0.001;
for epoch = 1:100
    xs(epoch+1) = xs(epoch) - alpha*2*xs(epoch);
end

figure;
plot(xs)