# vim:fdm=indent
import os
import os.path
import re
import glob
import subprocess
import sys
from shutil import copy2, copyfileobj, rmtree

def get_file_list(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError("file %s does not exist" % directory)

    return os.listdir(directory)

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

def retrieve(directory, file_list):
    if not os.path.isdir(directory):
        raise FileNotFoundError("file %s does not exist" % directory)

    for fl in file_list:
        copy2(directory + '/' + fl, './' + fl)

def copy(directory, file_list, new_directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError("directory %s does not exist" % directory)

    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)

    for fl in file_list:
        copy2(directory + '/' + fl, new_directory + '/' + fl)

def compress(tarfile, directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError("directory %s does not exist" % directory)

    returncode = subprocess.call(['tar', '-cf', tarfile, directory])

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
