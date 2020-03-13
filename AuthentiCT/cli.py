"""deam2cont -- Use deamination patterns to estimate contamination rates

"""

from argparse import ArgumentParser, FileType, RawTextHelpFormatter

import sys
import numpy as np
import numdifftools as ndt
import pandas as pd
import random
from scipy.optimize import minimize


from .input import *
from .contamination_est import compute_loglikelihood, format_loglikelihood
from .simulation import *

def probability(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def get_configuration():

    parent_parser = ArgumentParser(add_help=False)

    parent_parser.add_argument('-o', dest='outputfile', metavar='FILE', type=FileType('w'),
                        help="Output file\n"
                             "(default: console)")

    parent_parser.add_argument('-c', '--config', metavar='FILE', type=FileType('r'), nargs=1,
                        help="Input configuration file in the following format:\n"
                        "e          0.003170\n"
                        "rss        0.602679\n"
                        "lo         0.261435\n"
                        "lss        0.142886\n"
                        "lds        0.008414\n"
                        "rds        0.035792\n"
                        "contam     0.001000\n"
                        "o          0.637999\n"
                        "o2         0.477368")

    parser = ArgumentParser(description='Tools to explore deamination patterns on ancient DNA',
                            formatter_class=RawTextHelpFormatter)

    subparsers = parser.add_subparsers(dest='command', help='sub-command help')




    parser_inference = subparsers.add_parser('deam2cont', parents=[parent_parser],
                        help='Use deamination patterns to estimate contamination rate',
                        formatter_class=RawTextHelpFormatter)

    parser_inference.add_argument('-t', '--terminal', action='store_true',
                        help="Estimate contamination rate from terminal C-to-T substitutions only\n"
                              "(default: False)")

    parser_inference.add_argument('-m', '--mapq', metavar='MQ', type=int, default=0,
                        help="Mapping quality cutoff\n"
                             "(default: %(default)s)")

    parser_inference.add_argument('-l', '--minlength', metavar='INT', type=int, default=0,
                        help="Read length cutoff\n"
                             "(default: %(default)s)")

    parser_inference.add_argument('-b', '--bq', metavar='BQ', type=int, default=0,
                        help="Base quality cutoff\n"
                             "(default: %(default)s)")

    parser_inference.add_argument('inputfile', metavar='FILE', type=FileType('r'), nargs=1,
                        help = "Input BAM file"
                               "(use `-` for STDIN)")

    parser_inference.add_argument('-p', '--positions', metavar='FILE', type=FileType('r'), nargs=1,
                        help = "Provide here the positions that sequences should overlap")

    parser_inference.add_argument('-s', '--sample', metavar='INT', type=int, default=100000,
                        help="Maximum number of sequences used to fit the deamination model\n"
                             "(default: %(default)s)")

    parser_inference.add_argument('--decoding', action='store_true',
                        help="Print the posterior probabilities of each state, one line per position\n"
                              "(default: False)")


    parser_deamination = subparsers.add_parser('deamination',
                        help="Print deamination patterns",
                        formatter_class=RawTextHelpFormatter)

    parser_deamination.add_argument('inputfile', metavar='FILE', type=FileType('r'), nargs=1,
                        help = "Input BAM file"
                               "(use `-` for STDIN)")

    parser_deamination.add_argument('-o', dest='outputfile', metavar='FILE', type=FileType('w'),
                        help="Output file\n"
                             "(default: console)")

    parser_deamination.add_argument('-d', '--distance', metavar='INT', type=int, #default=5,
                        help="Compute the observed and expected distance between internal deaminations\n"
                             "Disregard the number of terminal positions provided\n"
                             "(default: %(default)s)")

    parser_deamination.add_argument('-m', '--mapq', metavar='MQ', type=int, default=0,
                        help="Mapping quality cutoff\n"
                             "(default: %(default)s)")

    parser_deamination.add_argument('-l', '--minlength', metavar='INT', type=int, default=0,
                        help="Read length cutoff\n"
                             "(default: %(default)s)")

    parser_deamination.add_argument('-b', '--bq', metavar='BQ', type=int, default=0,
                        help="Base quality cutoff\n"
                             "(default: %(default)s)")



    parser_simulation = subparsers.add_parser('simulation', parents=[parent_parser],
                        help='Simulate ancient DNA sequences',
                        formatter_class=RawTextHelpFormatter)

    parser_simulation.add_argument('-N', metavar='INT', type=int, default=10000,
                        help="Number of simulated sequences to generate\n"
                             "(default: %(default)s)")

    parser_simulation.add_argument('-GC', '--GCcontent', metavar='FLOAT', type=probability, default=0.4,
                        help="GC rate in simulated sequences\n"
                             "(default: %(default)s)")

    parser_simulation.add_argument('-L', '--Length', metavar='INT', type=int, default=50,
                        help="Average length of simulated sequences\n"
                             "(default: %(default)s)")

    parser_simulation.add_argument('-l', '--minlength', metavar='INT', type=int, default=0,
                        help="Minimum length of simulated sequences\n"
                             "(default: %(default)s)")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    opts = parser.parse_args()
    if opts.command == 'deam2cont' or opts.command == 'deamination':
        opts.inputfile = opts.inputfile[0]
    return opts


def main():

    options = get_configuration()

    if options.command == "simulation":
        if options.config:
            parameters = np.array(read_parameters(options.config[0]))
        else:
            parameters = np.array([0.002, 0.87, 0.35, 0.21, 0.003, 0.015, 0.001, 0.59, 0.49])
        
        for outputline in format_simulated_reads(parameters, options.Length, options.minlength, options.GCcontent, options.N):
            print(*outputline, sep='\t', file=options.outputfile)
        return

    elif options.command == "deamination":
        reads = list(read_input(options.mapq, options.minlength, options.bq, fh=options.inputfile))
        if options.distance:
            print("#distance\tobserved_rate\texpected_rate\tcount", file=options.outputfile)
            for outputline in distance_btw_deam(reads, options.distance):
                print(*outputline, sep='\t', file=options.outputfile)
        else:
            print("#position\trate\tderived\ttotal\tlower_95%CI\tupper_95%CI", file=options.outputfile)
            for outputline in deamination_patterns(reads):
                print(*outputline, sep='\t', file=options.outputfile)

    elif options.command == "deam2cont":
        if options.positions:
            all_sites = read_sites(options.positions[0])
            if options.config or options.terminal == True:
                reads = list(read_input(options.mapq, options.minlength, options.bq, all_sites, filter_reads=True, fh=options.inputfile))
            else:
                reads = list(read_input(options.mapq, options.minlength, options.bq, all_sites, fh=options.inputfile))
        else:
            reads = list(read_input(options.mapq, options.minlength, options.bq, fh=options.inputfile))


        if options.terminal == True:
            fragments = count_terminal_CtoT(reads)
            if 0 in fragments.values():
                print("Could not find DNA fragments in the following categories:")
                print(list(fragments.keys())[list(fragments.values()).index(0)])
                print("nodeam:", fragments["nodeam"], "; deam5:", fragments["deam5"], "; deam3:", fragments["deam3"], "; deam53:", fragments["deam53"], sep='\t')
                return
            total = fragments["nodeam"] + fragments["deam5"] + fragments["deam3"] + fragments["deam53"]
            contamination = fragments["nodeam"] / total - fragments["deam5"] * fragments["deam3"] / (fragments["deam53"] * total)
            print("#contamination_rate\t#nodeam\t#deam5\t#deam3\t#deam53", file=options.outputfile)
            print(contamination, fragments["nodeam"], fragments["deam5"], fragments["deam3"], fragments["deam53"], sep='\t', file=options.outputfile)
        else:
            if options.config:
                param = np.array(read_parameters(options.config[0]))
            else:
                guess = np.array([0.001, 0.9, 0.4, 0.003, 0.01, 0.01, 0.5, 0.55, 0.55])
                bnds = ((0.0001, 0.1),
                        (0.001, 0.999),
                        (0.001, 0.999),
                        (0.001, 0.999),
                        (0.001, 0.999),
                        (0.001, 0.1),
                        (0.001, 0.999),
                        (0.01, 0.99),
                        (0.01, 0.99))

                sample_size = options.sample
                sample = [reads[i] for i in random.sample(range(len(reads)), min(sample_size, len(reads)))]

                res = minimize(compute_loglikelihood, x0=guess, args=sample,
                           method='L-BFGS-B', bounds=bnds, tol=1e-10)

                Hfun = ndt.Hessian(compute_loglikelihood, step=0.0001, full_output=True)
                hessian_ndt, info = Hfun(res['x'], sample)
                with np.errstate(invalid='raise'):
                    try:
                        se = np.sqrt(np.diag(np.linalg.inv(hessian_ndt)))
                    except FloatingPointError:
                        print("It was not possible to compute the standard error for some estimates; perhaps you did not use enough sequences.")
                        with np.errstate(invalid='ignore'):
                            se = np.sqrt(np.diag(np.linalg.inv(hessian_ndt)))
                results = pd.DataFrame({'parameters':res['x'],'std err':se})
                results.index=['e','rss','lo','lss','lds','rds','contam','o','o2']

                param = res.x
                print(results, file=options.outputfile, flush=True)


            for outputline in format_loglikelihood(param, reads, per_position=options.decoding):
                print(*outputline, sep='\t', file=options.outputfile)
    else:
        print("Unrecognized command line argument")
if __name__ == '__main__':
    main()

