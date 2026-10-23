"""Entrypoint for the pinned modern-quantifier comparison."""
import argparse
import subprocess
import sys
from pathlib import Path
from benchmark_common import write_json
from .environment import bootstrap, load_environment, tool_paths, read_json
from .reference import canonical_reference, official_nac_sources
from .inputs import bulk, single_cell, single_end, prefix_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=Path('target/benchmarks/modern-comparison.noindex'))
    parser.add_argument('--data', type=Path, default=Path('data'))
    parser.add_argument('--tool', action='append', default=[], metavar='NAME=PATH')
    parser.add_argument('action', choices=['bootstrap', 'retrieve', 'prepare', 'plan', 'run', 'report', 'pilot', 'all'])
    parser.add_argument('--reference', choices=['gencode50', 'official-nac'], default='gencode50')
    parser.add_argument('--single-cell-strand', choices=['fw', 'both'], default='fw')
    parser.add_argument('--amend-unstarted-workflow', action='store_true',
                        help='Audit a workflow revision after the index panel, before any read measurements')
    parser.add_argument('--only', help='Comma-separated case IDs, or index')
    args = parser.parse_args()
    root, data = args.workspace.resolve(), args.data.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if args.action == 'bootstrap':
        bootstrap(root, dict(item.split('=', 1) for item in args.tool))
    elif args.action == 'all':
        bootstrap(root, dict(item.split('=', 1) for item in args.tool))
        # Build once, validate commands on prefixes, then collect the read panel.
        actions = [('retrieve', None), ('prepare', None), ('run', 'index')]
        if args.only != 'index':
            read_selection = ','.join(v for v in args.only.split(',') if v != 'index') if args.only else None
            actions += [('pilot', read_selection), ('run', read_selection)]
        actions += [('report', None)]
        for action, selection in actions:
            command = [sys.executable, str(Path(__file__).resolve().parents[1] / 'modern_benchmarks.py'),
                       '--workspace', str(root), '--data', str(data), '--reference', args.reference, '--single-cell-strand', args.single_cell_strand, action]
            if selection:
                command += ['--only', selection]
            subprocess.run(command, check=True)
    elif args.action == 'retrieve':
        env = tool_paths(load_environment(root))
        scripts = Path(__file__).resolve().parents[1]
        for accession in (('SRR10919051',) if args.reference == 'official-nac' else ('SRR13638690', 'SRR307903')):
            if (data / accession.lower() / 'trimmed-fastp/preparation.json').exists():
                bulk(data, accession)
            else:
                subprocess.run([sys.executable, str(scripts / 'prepare_benchmark.py'), '--data', str(data),
                                '--fastp', env['fastp'], '--threads', '8', '--samples', accession], check=True)
        if args.reference == 'official-nac':
            official_nac_sources(data)
        else:
            if (data / 'tenx-pbmc-1k-v3/preparation.json').exists():
                single_cell(data)
            else:
                subprocess.run([sys.executable, str(scripts / 'prepare_benchmark.py'), '--data', str(data),
                                '--fastp', env['fastp'], '--threads', '8', '--tenx-pbmc-v3'], check=True)
            from benchmark_common import fetch
            fetch('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz',
                  data / 'srr-small/gencode.v50.transcripts.fa.gz',
                  sha256='5a320f524d73b5793518eb19b118829033713443d0f42af20a67bb31cc06cf56', size=183554921)
    elif args.action == 'prepare':
        env = load_environment(root)
        if args.reference == 'official-nac':
            reference = canonical_reference(official_nac_sources(data), root / 'reference', kind='official-nac')
            cases = [bulk(data, 'SRR10919051')]
            cases += [single_end(case) for case in cases.copy()]
        else:
            reference = canonical_reference(data / 'srr-small/gencode.v50.transcripts.fa.gz', root / 'reference')
            cases = [bulk(data, x) for x in ('SRR13638690', 'SRR307903')]
            from .simulation import simulate
            cases.append(simulate(reference, root / 'truth'))
            cases += [single_end(case) for case in cases.copy()]
            sc = single_cell(data)
            sc['strand'] = args.single_cell_strand
            cases.append(sc)
        write_json(root / 'inputs.json', {'reference': reference, 'cases': cases})
        print('Prepared and verified reference and libraries', flush=True)
    elif args.action in ('run', 'plan', 'pilot'):
        from .experiment import execute, plan
        env = load_environment(root)
        inputs = read_json(root / 'inputs.json')
        if args.action == 'plan':
            print(plan(inputs))
        elif args.action == 'pilot':
            pilot = {**inputs, 'cases': [prefix_case(c, root / 'pilot-inputs' / c['id']) for c in inputs['cases']]}
            selection = ','.join(v if v == 'index' or v.endswith('-pilot') else v + '-pilot'
                                 for v in args.only.split(',')) if args.only else None
            execute(root, tool_paths(env), pilot, pilot=True, only=selection)
        else:
            execute(root, tool_paths(env), inputs, only=args.only, amend=args.amend_unstarted_workflow)
    else:
        from .tables import report
        report(root)


if __name__ == '__main__':
    main()
