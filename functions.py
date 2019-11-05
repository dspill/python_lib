# vim:fdm=indent
import os
import os.path
import re
import glob
import subprocess
import contextlib
import io
import sys
import numpy as np
from shutil import copy2, copyfileobj

# Python3
def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                func(*a, **ka)
    return wrapper

# Python2
@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout

def YesNo(Question):
    """Ask for yes or no answer and return a boolean."""

    Yes = set(["YES", "Y", "yes", "y", "Yes", ""])
    No = set(["NO", "N", "no", "n", "No"])

    while True:
        Answer = input(Question).lower()
        if Answer in Yes:
            return True
        elif Answer in No:
            return False
        else:
            print("Possible answers:")
            print("  %s" % sorted(list(Yes)))
            print("or")
            print("  %s" % sorted(list(No)))


def read_xy(datafile, xcol=0, ycol=1):
    return np.loadtxt(datafile, unpack=True, usecols=(xcol, ycol))

''' write multi-column file '''
def write_multiple(list_of_lists, filename):
    ll = np.array(list_of_lists).transpose()

    with open(filename, 'w') as outfile:
        for row in ll:
            for column in row:
                outfile.write('%15.9e ' % column)
            outfile.write('\n')

    outfile.close()

def write_xy(x, y, filename):
    if len(x) != len(y):
        raise RuntimeError("lengths of lists is not compatible")

    write_multiple([x, y], filename)



def compactifyDump(filename="./dump.tar", compress=False):
    # find largest step
    files = glob.glob('./dump/fluid*.0*')
    if len(files) == 0:
        start_step = 0
    else:
        steps = []
        for file in files:
            new = re.search('.+fluid(.+)\.0\.dat', file)
            steps.append(int(new.group(1)))

        steps.sort()
        start_step = max(steps)

    print("largest step = " + str(start_step))

    # copy most recent output to separate directory
    if not os.path.isdir("./dump_last"):
        print("creating dump_last")
        os.mkdir("./dump_last")

    files = glob.glob("./dump/*"+str(start_step)+".*")
    print("copying files")
    for file in files:
        copy2(file, "./dump_last/")

    # compress remaining output and remove directory if successful
    if os.path.isdir("./dump"):
        print("compressing dump")
        if compress:
            if subprocess.call(['tar', '-cf', filename, '-I', 'pigz', './dump']) == 0:
                print("compression successful, removing dump")
                # rmtree("./dump")
        else:
            if subprocess.call(['tar', '-cf', filename, './dump']) == 0:
                print("taring successful, removing dump")
                # rmtree("./dump")

def concatenateFiles(outfile, file_list):
    with open(outfile, 'wb') as wfd:
        for f in file_list:
            with open(f, 'rb') as fd:
                copyfileobj(fd, wfd, 1024*1024*10)
                #10MB per writing chunk to avoid reading big file into memory.


def get_file_list_from_tar(tarfile):
    if not os.path.isfile(tarfile):
        raise FileNotFoundError("file %s does not exist" % tarfile)

    output = subprocess.check_output(['tar', '-tf', tarfile])
    output = output.decode("utf-8").split('\n')
    return output


def get_steps(file_list):
    steps = []
    for filename in file_list:
        z = re.match('.*fluid(.+)\.0\.dat', filename)
        if z:
            steps.append(int(z.groups()[0]))

    steps.sort()
    return steps


def reduce_file_list(file_list, name, step):
    new_list = []
    for filename in file_list:
        z = re.match('.*' + name + str(step) + '\.\d+\.dat', filename)

        if z:
            new_list.append(filename)

    new_list.sort()
    return new_list


def extract(tarfile, file_list):
    if not os.path.isfile(tarfile):
        raise FileNotFoundError("file %s does not exist" % tarfile)

    returncode = subprocess.call(['tar', '-xvf', tarfile] + file_list)

    if returncode != 0:
        raise ChildProcessError("tar command failed")


def concatenateFiles(file_list, outfile):
    with open(outfile, 'wb') as wfd:
        for f in file_list:
            with open(f, 'rb') as fd:
                #10MB per writing chunk to avoid reading big file into memory
                copyfileobj(fd, wfd, 1024*1024*10)


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def hours_minutes_seconds(td):
    ''' works for datetime.timedelta object '''
    return '{:02}:{:02}:{:02}'.format(
            int(td.days * 24 + td.seconds//3600),
            int(td.seconds//60)%60,
            int(td.seconds%60))
