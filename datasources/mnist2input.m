function [inp, ts, lbls] = mnist2input(seconds)

[patterns, labels] = loadmnistasspikes();
res = [];
lbls = [];
offset = 0;
seperation = 500;
steps = ceil(seconds * 1000 / seperation);
perm = randperm(10000, steps);

for i = 1 : steps
    patt = patterns{perm(i)};
    lbl = labels(perm(i));
    patt(:, 2) = patt(:, 2) + offset;
    res = [res; patt];
    lbls = [lbls; lbl];
    offset = offset + seperation;
end

[inp, ts] = sortspiketimes(res(:, 1)', res(:, 2)');

end