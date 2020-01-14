exp = {'jit', 'fgi', 'frq', 'naf', 'drp'};
expname = {'Internal pattern jitter', 'FGI', 'Pattern presentation frequency', 'Number of afferents', 'Pattern dropout'};
xlabels = {'Jitter (ms)', 'FGI', 'Frequency (Hz)', '# Afferents', 'Dropout percentage'};
rule = {'ss', 'sd'};
template = 'results/latest/res_%s%s_final.mat';

for i = 1 : numel(exp)
    for j = 1 : numel(rule)
        filename = sprintf(template, exp{i}, rule{j});
        load(filename);
        plotres;
        
        xlabel(xlabels{i});
        rulename = 'sSDVL';
        if strcmp(rule{j}, 'sd')
            rulename = 'SDVL';
        end
        title(sprintf('%s %s',  expname{i}, rulename));
        savename = sprintf('%stpxtn%s.eps', rulename, exp{i})
        saveas(figh, savename, 'epsc');
    end
end