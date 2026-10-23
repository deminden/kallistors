"""Explicit native command plans. No quantifier implementation is modified."""
from pathlib import Path

BULK = ('kallistors', 'salmon-sketch', 'salmon-sa', 'piscem-infer')
SC = ('kallistors-bustools', 'piscem-alevin-fry')
BUILDERS = ('kallistors', 'salmon', 'piscem')


def index_path(root, method):
    return root / ('salmon' if method.startswith('salmon') else
                   'piscem/index' if method.startswith('piscem') else 'kallistors.idx')


def index_plan(tools, builder, reference, root, threads):
    fasta = reference['fasta']['path']
    out = index_path(root, builder)
    if builder == 'kallistors':
        cmd = [tools[builder], 'index', '-i', str(out), '-t', str(threads), fasta]
    elif builder == 'salmon':
        cmd = [tools[builder], 'index', '-t', fasta, '-i', str(out), '-p', str(threads),
               '-k', '31', '-m', '19', '--keepDuplicates', '--no-clip']
    else:
        cmd = [tools[builder], 'build', '-s', fasta, '-o', str(out), '-t', str(threads),
               '-k', '31', '-m', '19', '--seed', '1', '--ram-limit-gib', '8',
               '--work-dir', str(root / 'scratch-dbg'), '--tmp-dir', str(root / 'scratch-sshash')]
    return [('index', threads, cmd)]


def bulk_plan(tools, method, case, reference, indexes, out, threads):
    idx = str(index_path(indexes, method))
    r1, r2 = ','.join(case['read1']), ','.join(case['read2'])
    t = str(threads)
    single = case.get('layout') == 'single'
    model = case.get('fragment_model', {})
    strand = case.get('strand', 'unstranded')
    library_type = ({'reverse': 'SR', 'forward': 'SF', 'unstranded': 'U'} if single else
                    {'reverse': 'ISR', 'forward': 'ISF', 'unstranded': 'IU'})[strand]
    if method == 'kallistors':
        cmd = [tools['kallistors'], 'quant', '-i', idx, '-o', str(out / 'quant'), '-t', t]
        if strand != 'unstranded':
            cmd.append('--rf-stranded' if strand == 'reverse' else '--fr-stranded')
        cmd += ['--single', '-l', str(model['mean']), '-s', str(model['sd']), *case['read1']] if single else sum(([a, b] for a, b in zip(case['read1'], case['read2'])), [])
        return [('quant', threads, cmd)]
    if method.startswith('salmon'):
        cmd = [tools['salmon'], 'quant', '-i', idx, '-o', str(out / 'quant'), '-p', t,
               '-l', library_type, '--sigDigits', '15']
        cmd += ['--fldMean', str(model['mean']), '--fldSD', str(model['sd']), '-r', *case['read1']] if single else ['-1', *case['read1'], '-2', *case['read2']]
        if method == 'salmon-sketch':
            cmd.append('--sketch')
        return [('quant', threads, cmd)]
    mapping = [tools['piscem'], 'map-bulk', '-i', idx, '-o', str(out / 'mapping'), '-t', t]
    mapping += ['-r', r1] if single else ['-1', r1, '-2', r2]
    inference = [tools['piscem-infer'], 'quant', '-i', str(out / 'mapping'),
                 '-o', str(out / 'quant'), '-l', library_type, '--num-threads', t]
    if single:
        inference += ['--fld-mean', str(model['mean']), '--fld-sd', str(model['sd'])]
    return [('mapping', threads, mapping), ('inference', threads, inference)]


def sc_plan(tools, method, case, reference, indexes, out, threads, mapping_only=False):
    idx, t = str(index_path(indexes, method)), str(threads)
    t2g = reference['gene_map']['path']
    mapped = out / 'mapping'
    if method == 'piscem-alevin-fry':
        plan = [('mapping', threads, [tools['piscem'], 'map-sc', '-i', idx, '-o', str(mapped),
                 '-t', t, '-g', 'chromium_v3', '-1', ','.join(case['read1']), '-2', ','.join(case['read2'])])]
        if threads == 1 or mapping_only:
            return plan
        permit, quant = out / 'permit', out / 'quant'
        return plan + [
            ('barcode_correction', threads, [tools['alevin-fry'], 'generate-permit-list', '-i', str(mapped),
             '-o', str(permit), '-d', case.get('strand', 'fw'), '-u', case['onlist'], '-m', '10', '-t', t,
             '--cell-bc-correction', 'unique', '--cell-bc-neighborhood', 'hamming-1']),
            ('collation', threads, [tools['alevin-fry'], 'collate', '-i', str(permit), '-r', str(mapped), '-t', t]),
            ('umi_resolution', threads, [tools['alevin-fry'], 'quant', '-i', str(permit), '-m', t2g,
             '-o', str(quant), '-t', t, '-r', 'cr-like'])]
    reads = sum(([a, b] for a, b in zip(case['read1'], case['read2'])), [])
    plan = [('mapping', threads, [tools['kallistors'], 'bus', '-i', idx, '-o', str(mapped),
                                  '-t', t, '-x', '10xv3', '--unstranded' if case.get('strand', 'fw') == 'both' else '--fr-stranded', *reads])]
    if mapping_only:
        return plan
    bt = tools['bustools']
    return plan + [
        ('sort_raw', threads, [bt, 'sort', '-t', t, '-m', '2G', '-o', str(out / 'raw.bus'), str(mapped / 'output.bus')]),
        ('barcode_correction', 1, [bt, 'correct', '-w', case['onlist'], '-o', str(out / 'corrected-unsorted.bus'), str(out / 'raw.bus')]),
        ('sort_corrected', threads, [bt, 'sort', '-t', t, '-m', '2G', '-o', str(out / 'corrected.bus'), str(out / 'corrected-unsorted.bus')]),
        ('umi_resolution', 1, [bt, 'count', '--genecounts', '-g', t2g, '-e', str(mapped / 'matrix.ec'),
         '-t', str(mapped / 'transcripts.txt'), '-o', str(out / 'quant'), str(out / 'corrected.bus')])]


def bulk_mapping_plan(tools, method, case, reference, indexes, out, threads):
    if method == 'piscem-infer':
        return None  # Its full pipeline already contains an isolated native mapping command.
    if method.startswith('salmon'):
        command = bulk_plan(tools, method, case, reference, indexes, out, threads)[0][2]
        return [('mapping', threads, command + ['--skipQuant'])]
    cmd = [tools['kallistors'], 'pseudoalign', '--index', str(index_path(indexes, method)),
           '--threads', str(threads), '--reads', case['read1'][0]]
    if case.get('strand', 'unstranded') != 'unstranded':
        cmd.append('--rf-stranded' if case['strand'] == 'reverse' else '--fr-stranded')
    if case.get('layout') == 'single':
        cmd += ['--fragment-length', str(case['fragment_model']['mean'])]
    else:
        cmd += ['--reads2', case['read2'][0]]
    return [('mapping', threads, cmd + ['--out', str(out / 'mapping.ec.tsv')])]
