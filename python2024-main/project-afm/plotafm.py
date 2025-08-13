from collections import Counter
from itertools import islice
from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt


def extract_data_points(fname):
    # return a dict called 'data', such that
    # data[s, i, j] = (d, f),
    # where s is the series, i, j are the point coordinates;
    # d is a numpy array containing measured distances,
    # f is a numpy array containing measured forces.
    
    data = dict()
    with open(fname, 'rt') as ftext:
        s, i, j = 0, 0, 0
        d, f = [], []
        reading_data = False

        for line in ftext:
            if line.startswith("#"):
                if "# index" in line:
                    if d and f and i is not None and j is not None:
                        data[(s, i, j)] = (np.array(d), np.array(f))
                        d, f = [], []
                        s = (s + 1) % 2 
                    reading_data = False
                elif "# iIndex" in line:
                    i = int(line.split()[-1])
                elif "# jIndex" in line:
                    j = int(line.split()[-1])
                elif "# units" in line:
                    reading_data = True
                continue

            if reading_data:
                values = line.split()
                if len(values) >= 2:
                    d.append(float(values[0]))
                    f.append(float(values[1]))

        if d and f and i is not None and j is not None:
            data[(s, i, j)] = (np.array(d), np.array(f))

    return data

def raw_plot(point, curve, save=None, show=True):
    """Plot one raw distance-force curve."""
    d, f = curve
    s, i, j = point
    plt.figure(figsize=[9, 6])
    plt.plot(d, f, label=f'Series {s}, Point ({i},{j})')
    plt.xlabel('Distance (m)')
    plt.ylabel('Force (N)')
    plt.title(f'Distance vs. Force for Series {s}, Point ({i},{j})')
    plt.legend()
    plt.grid(True)
    if save is not None:
        plt.savefig(save, dpi=200, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()

def do_raw_plots(data, show, plotprefix):
    for point, curve in data.items():
        s, i, j = point
        print(f"Plotting curve at {point}")
        fname = f'{plotprefix}-{s:01d}-{i:03d}-{j:03d}.png' if plotprefix is not None else None
        raw_plot(point, curve, show=show, save=fname)

def main(args):
    fname = args.textfile
    print(f"Parsing {fname}...")
    full_data = extract_data_points(fname)
    if args.first is not None:
        data = dict((k, v) for k, v in islice(full_data.items(), args.first))
    else:
        data = full_data
    do_raw_plots(data, args.show, args.plotprefix)

def get_argument_parser():
    p = ArgumentParser()
    p.add_argument("--textfile", "-t", required=True,
                   help="Name of the data file containing AFM curves for many points")
    p.add_argument("--first", type=int,
                   help="Number of curves to extract and plot")
    p.add_argument("--plotprefix", default="curve",
                   help="Non-empty path prefix of plot files (PNGs); do not save plots if not given")
    p.add_argument("--show", action="store_true",
                   help="Show each plot")
    return p

if __name__ == "__main__":
    p = get_argument_parser()
    args = p.parse_args()
    main(args)

