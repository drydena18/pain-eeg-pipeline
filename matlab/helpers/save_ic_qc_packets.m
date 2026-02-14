function save_ic_qc_packets(QC, subjid, EEG, icList)
if isempty(icList), return; end
ensure_dir(QC);

outDir = fullfile(QC, sprintf('sub-%03d_icqc', subjid));
ensure_dir(outDir);

fs = EEG.srate;

W = EEG.icaweights * EEG.icasphere;
X = double(EEG.data(EEG.icachansind, :));
act = W * X;

C = [];
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classifications') && isfield(EEG.etc.ic_classifications, 'ICLabel') && isfield(EEG.etc.ic_classifications.ICLabel, 'classifications')
    C = EEG.etc.ic_classification.ICLabel.classifications;
end

for ic = icList(:)'
    h = figure('Visible', 'off');

    subplot(2, 2, 1);
    if isfield(EEG, 'icawinv') && has_chanlocs(EEG)
        topoplot(EEG.icawinv(:, ic), EEG.chanlocs, 'electrodes', 'on');
        title(sprintf('IC %d scalp map', ic), 'Interpreter', 'none');
    end

    subplot(2, 2, 2);
    win = round(fs * 2);
    nover = round(win * 0.5);
    nfft = max(2^nextpow2(win), win);
    [pxx, f] = pwelch(act(ic, :), win, nover, nfft, fs);
    plot(f, 10*log10(pxx));
    xlim([0 80]);
    xlabel('Hz'); ylabel('Power (dB)');
    title('IC PSD', 'Interpreter', 'none');

    subplot(2, 2, [3 4]);
    nSamp = min(size(act, 2), round(fs * 10));
    plot((0:nSamp-1) / fs, act(ic, 1:nSamp));
    xlabel('Time (s)'); ylabel('a.u.');
    title('IC activation (first 10s)', 'Interpreter', 'none');

    if ~isempty(C) && size(C, 1) >= ic
        p = C(ic, :);
        ttl = sprintf('sub-%03d IC %d | B%.2f M%.2f E%.2f H%.2f L%.2f C%.2f O%.2f', ...
            subjid, ic, p(1), p(2), p(3), p(4), p(5), p(6), p(7));
    else
        ttl = sprintf('sub-%03d IC %d', subjid, ic);
    end
    sgtitle(ttl, 'Interpreter', 'none');

    saveas(h, fullfile(outDir, sprintf('sub-$03d_ic%03d_qc.png', subjid, ic)));
    close(h);
end
end