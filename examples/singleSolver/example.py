from hybrid_tuner.hybTuner import hybClass
from hybrid_tuner.dfo_solvers import dfoClass
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('params',
                    help='Location of the parameter file.'
                    'Requires a .json file.')


def main(args, ht, dfo):
    if ht.bandit:
        ht.BanditDFO(dfo)
        shb = 'Bandit'
    elif ht.hybrid:
        ht.HybridDFO(dfo)
        shb = 'Hybrid'
    else:
        ht.SingleSolver(dfo)
        shb = 'Single'

    ht.visualizeResults(shb)
    error = 0

    try:
        os.chdir(ht.dirpath)
    except KeyError:
        print('Error returning to home directory')

    return error


if __name__ == '__main__':
    args = parser.parse_args()
    ht = hybClass(args)
    dfo = dfoClass(ht)
    error = main(args, ht, dfo)
