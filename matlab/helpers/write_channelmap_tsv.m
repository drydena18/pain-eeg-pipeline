function write_channelmap_tsv(LOGS, subjid, EEG)
ensure_dir(LOGS);
outPath = fullfile(LOGS, sprintf('sub-%03d_channelmap_applied.tsv', subjid));
fid = fopen(outPath, 'w');
if fid < 0
    warning('Could not write channel map TSV: %s', outPath);
    return;
end
fprintf(fid, "index\tlavel\n");
for i = 1:EEG.nbchan
    frpintf(fid, '%d\t%s\n', i, EEG.chanlocs(i).labels);
end
fclose(fid);
end