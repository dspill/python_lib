import argparse
import re
import os
import time
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

 # value given by Latex \the\columnwidth
columnwidths = {
        'beamer' : 307.28987
        }

def filter_files(data_dir, regexp):
    # regexp = '.+thr-1.+'
    regexp_comp = re.compile(regexp)

    filenames = sorted(os.listdir(data_dir))
    filenames = np.array(list(filter(regexp_comp.match, filenames)))
    lst = [data_dir + '/' + filename for filename in filenames]
    return sorted(lst)


def read_xy(datafile, xcol=0, ycol=1):
    return np.loadtxt(datafile, unpack=True, usecols=(xcol, ycol))

def read_cols(datafile, columns):
    return np.loadtxt(datafile, unpack=True, usecols=columns)

def read_xyz(filename):
    with open(filename) as f:
        line = f.readline()
        num_particles = int(line.split()[0])
        f.readline() # skip box

        pid  = np.zeros(num_particles, dtype=int)
        xpos = np.zeros(num_particles)
        ypos = np.zeros(num_particles)
        zpos = np.zeros(num_particles)

        for i in range(num_particles):
            line = f.readline().split()

            pid[i]  = int(line[0])
            xpos[i] = float(line[1])
            ypos[i] = float(line[2])
            zpos[i] = float(line[3])

    return pid, xpos, ypos, zpos


def YesNo(Question):
    """Ask for yes or no answer and return a boolean."""

    Yes = set(["YES", "Y", "yes", "y", "Yes", ""])
    No = set(["NO", "N", "no", "n", "No"])

    while True:
        Answer = input(Question).lower()
        if Answer in Yes:
            return True
        if Answer in No:
            return False

        print("Possible answers:")
        print("  %s" % sorted(list(Yes)))
        print("or")
        print("  %s" % sorted(list(No)))


def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()


def trim(x, y, xmin, xmax):
    assert len(x) == len(y)
    xn = []
    yn = []
    for i, val in enumerate(x):
        if xmin <= val <= xmax:
            xn.append(val)
            yn.append(y[i])
    return xn, yn

def parabola(x, a, x0, y0):
    return a * (x - x0)**2 + y0

def fit_peak0(x, y, points):
    assert len(x) >= 3
    points = sorted(points, key=lambda l: l[0])

    x0 = points[1][0]
    y0 = points[1][1]

    a = (points[0][0] - x0) / (points[0][1] - y0)**2
    a += (points[2][0] - x0) / (points[2][1] - y0)**2
    a /= 2.

    p0 = [a, x0, y0]
    popt, pcov = curve_fit(parabola, x, y, p0=p0)

    return popt, pcov


def fit_peak(x, y, axes, filename, step):
    tellme('Select three points to draw a triangle around the peak')

    happy = False
    while not happy:
        pts = []
        while len(pts) < 3:
            tellme('Select 3 corners with mouse')
            pts = np.asarray(plt.ginput(3, timeout=-1))
            if len(pts) < 3:
                tellme('Too few points, starting over')
                time.sleep(1)  # Wait a second

        xmin = min(pts[:,0])
        xmax = max(pts[:,0])
        xn, yn = trim(x, y, xmin, xmax)

        try:
            popt, pcov = fit_peak0(xn, yn, pts)
            print("x_p = " + str(popt[1]))
            xls = np.linspace(xmin, xmax, 200)
            axes.plot(xls, parabola(xls, *popt))
            plt.draw()
            tellme('Happy? Key press for yes, mouse click for no')
            happy = plt.waitforbuttonpress()
        except RuntimeError:
            print('Could not find optimal parameters')

        if happy:
            # write to file
            outfile_name = filename
            if os.path.isfile(outfile_name):
                print('File ' + outfile_name
                        + ' already exists. Appending step '
                        + str(step))

            outfile = open(outfile_name, 'a')
            outfile.write('%10d %16.9e %16.9e\n'
                    % (step, popt[1], pcov[1][1]))
            outfile.close()
            plt.cla()
        else:
            tellme('Abort? Key press for yes, mouse click for no')
            if plt.waitforbuttonpress():
                plt.cla()
                break


def select_point(x, y, axes, outfile_name, step):
    tellme('Select point')
    pt = plt.ginput()[0]
    print('selected ', pt)
    # write to file
    if os.path.isfile(outfile_name):
        print('File ' + outfile_name
                + ' already exists. Appending step '
                + str(step))

    outfile = open(outfile_name, 'a')
    outfile.write('%10d %16.9e %16.9e\n' % (step, pt[0], pt[1]))
    outfile.close()
    plt.cla()


def select_point_xerr(x, y, axes, outfile_name, step):
    tellme('Select point')
    pts = np.asarray(plt.ginput(3, timeout=-1))
    pts = sorted(pts , key=lambda k: k[0]) # sort according to x (0) coordinate
    print('selected ', pts)
    assert(pts[0][0] <= pts[1][0] <= pts[2][0])

    # write to file
    if os.path.isfile(outfile_name):
        print('File ' + outfile_name
                + ' already exists. Appending step '
                + str(step))

    outfile = open(outfile_name, 'a')
    outfile.write('%10d %16.9e %16.9e %16.9e %16.9e\n' % (step, pts[1][0],
        pts[1][1], pts[0][0], pts[2][0]))
    outfile.close()
    plt.cla()


def select_point_yerr(x, y, axes, outfile_name, step):
    tellme('Select point')
    pts = np.asarray(plt.ginput(3, timeout=-1))
    pts = sorted(pts , key=lambda k: k[1]) # sort according to y (1) coordinate
    print('selected ', pts)
    assert(pts[0][1] <= pts[1][1] <= pts[2][1])

    # write to file
    if os.path.isfile(outfile_name):
        print('File ' + outfile_name
                + ' already exists. Appending step '
                + str(step))

    outfile = open(outfile_name, 'a')
    outfile.write('%10d %16.9e %16.9e %16.9e %16.9e\n' % (step, pts[1][0],
        pts[1][1], pts[0][0], pts[2][0]))
    outfile.close()
    plt.cla()


def ParseArguments():
    parser = argparse.ArgumentParser(description='''plot your data files''')

    parser.add_argument('-f', '--input', nargs='+', type=str,
            help='list of input files', required=True)

    parser.add_argument('-r', '--from_row', action='store_true',
            help='Take row index as x-axis')

    parser.add_argument('-m', '--marker', action='store_true',
            help='use markers in plot')

    parser.add_argument('-x', '--xcol', type=int, default=0,
            help='x-column to plot')

    parser.add_argument('-y', '--ycol', nargs='+', type=int,
            default=[1], help='y-column(s) to plot')

    parser.add_argument('--xfac', type=int, default=1,
            help='factor by which to scale x-axis')

    parser.add_argument('-o', '--output', type=str, help='output file')

    parser.add_argument('--fit_peak', type=str,
            help='fit parabola to peak and store position in given file')

    parser.add_argument('--select_point', type=str,
            help='select point in plot and store it in given file')

    parser.add_argument('--select_point_xerr', type=str,
            help='select three points in plot, which form the actual data '
            'point plus error bars and store it in given file')

    parser.add_argument('--select_point_yerr', type=str,
            help='select three points in plot, which form the actual data '
            'point plus error bars and store it in given file')

    parser.add_argument('--slice', action='store_true',
            help='plot a sliced (2d) configuration')

    parser.add_argument('--logx', action='store_true', help='logarithmic x-axis')
    parser.add_argument('--logy', action='store_true', help='logarithmic y-axis')

    parser.add_argument('--xmax', type=float, help='upper limit x-axis')
    parser.add_argument('--ymax', type=float, help='upper limit y-axis')

    parser.add_argument('--xmin', type=float, help='lower limit x-axis')
    parser.add_argument('--ymin', type=float, help='lower limit y-axis')

    return parser.parse_args()


def get_step(string):
    z = re.match('\D*0*(\d+)\D*', string)
    if z:
        # print('step = ' + str(int(z.groups()[0])))
        return int(z.groups()[0])

    return int(input('timestep? '))


def get_integer(string):
    z = re.match('\D*(\d+)\D*', string)
    if z:
        return int(z.groups()[0])

    print("no integer detected for " + string)
    return None


def get_decimal(string):
    z = re.match('.*(\d+[.]\d+)\D*', string)
    if z:
        return float(z.groups()[0])


def plot_slice(datafiles):
    n_files = len(datafiles)
    n_cols = int(round(math.sqrt(n_files)))

    # do the plotting
    _, axes = plt.subplots(math.ceil(n_files / n_cols), n_cols)

    if n_files == 1:
        lst = [axes]
    else:
        lst = axes.flatten()

    for i, axis in enumerate(lst):
        axis.set_aspect('equal')
        if i < n_files:
            (_, xpos, ypos, _) = read_xyz(datafiles[i])
            axis.scatter(xpos, ypos, marker='o', s=2, linewidth=0., label=datafiles[i])

            title = datafiles[i].replace("sliced_frame_", "t = ")
            title = title.replace(".xyz", "")
            axis.set_title(title, fontsize=6)
            # axis.legend(loc=1)

        # disable ticks
        axis.tick_params(
                axis='both',       # changes apply to both axes
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                left=False,
                labelbottom=False,
                labelleft=False,
                )

        plt.tight_layout()
    # plt.axes().set_aspect('equal', 'datalim')
    plt.savefig('configurations.pdf')
    plt.show()

def plot(args):
    datafiles = args['input']
    # do the plotting
    if args['slice']:
        plot_slice(datafiles)
        return

    _, axes = plt.subplots()

    if args['logx']:
        axes.set_xscale('log')
    if args['logy']:
        axes.set_yscale('log')

    for i_f, f in enumerate(datafiles):
        for ycol in args['ycol']:
            x, y = read_xy(f, args['xcol'], ycol)

            if args['from_row']:
                x = [i for i in range(len(y))]

            if args['marker']:
                marker = '.'
            else:
                marker = None

            x *= args['xfac']
            axes.plot(x, y, label=f+' c'+str(ycol), marker=marker)

        if i_f == 0:
            min_x = min(x)
            max_x = max(x)
            min_y = min(y)
            max_y = max(y)
        else:
            min_x = min(min_x, min(x))
            max_x = max(max_x, max(x))
            min_y = min(min_y, min(y))
            max_y = max(max_y, max(y))

        dx = .05*(max_x - min_x)
        dy = .05*(max_y - min_y)
        xlim = np.array([min_x - dx, max_x + dx])
        ylim = np.array([min_y - dy, max_y + dy])

        if args['xmin']:
            xlim[0] = args['xmin']
        if args['xmax']:
            xlim[1] = args['xmax']

        if args['ymin']:
            ylim[0] = args['ymin']
        if args['ymax']:
            ylim[1] = args['ymax']

        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

        axes.legend()

        if args['fit_peak']:
            fit_peak(x, y, axes, args['fit_peak'], get_step(f))

        if args['select_point']:
            select_point(x, y, axes, args['select_point'], get_step(f))

        if args['select_point_xerr']:
            select_point_xerr(x, y, axes, args['select_point_xerr'], get_step(f))
        if args['select_point_yerr']:
            select_point_yerr(x, y, axes, args['select_point_yerr'], get_step(f))


    if args['output']:
        plt.savefig(args['output'])
    if args['fit_peak']:
        plt.close()
    else:
        plt.show()


def get_figsize(wf=0.5, hf=(5.**0.5-1.0)/2.0, style='beamer'):
    """Parameters:
      - wf [float]:     width fraction in columnwidth units
      - hf [float]:     height fraction in columnwidth units.
                        Set by default to golden ratio.
      - style [string]: Type of document. Get columnwidth from LaTeX
                        using \showthe\columnwidth
    Returns:  (fig_width,fig_height): that should be given to matplotlib
    """
    if style is 'beamer': 
        columnwidth = 307.28987 # pt
    elif style is 'article':
        columnwidth = 455.24411 # pt
    else:
        raise ValueError('style %s not defined' % style)

    fig_width_pt = columnwidth*wf
    inches_per_pt = 1.0/72                 # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*hf      # height in inches
    return (fig_width, fig_height)


def toggle_selector(event):
    print('Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print('RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print('RectangleSelector activated.')
        toggle_selector.RS.set_active(True)
