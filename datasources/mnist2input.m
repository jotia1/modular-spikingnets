function [inp, ts] = mnist2input(seconds)

patterns = loadmnistasspikes();
res = [];
offset = 0;
seperation = 500;
steps = ceil(seconds * 1000 / seperation);
perm = randperm(10000, steps);

for i = 1 : steps
    patt = patterns{perm(i)};
    patt(:, 2) = patt(:, 2) + offset;
    res = [res; patt];
    offset = offset + seperation;
end

plot(res(:, 2), res(:, 1), '.k'); 
ax = gca;
ax.YDir = 'normal';










[inp, ts] = sortspiketimes(res(:, 1)', res(:, 2)');

end