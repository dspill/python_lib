#!/usr/bin/env python3
import os
import re
import tarfile
import tempfile
import subprocess
import warnings
import shutil
from shutil import copy2, copyfileobj
from distutils.dir_util import copy_tree

class Dump_container:
    def __init__(self, file_list):
        self.names = ['pops', 'fluid', 'couplForces']
        self.file_list = file_list

        # get list of steps
        self.steps = []
        for filename in file_list:
            z = re.match('.*fluid(.+)\.0\.dat', filename)
            if z:
                self.steps.append(int(z.groups()[0]))

        self.steps.sort()
        self.check_consistency()

    def n_cores(self):
        return len(self.reduced_file_list(steps=(self.steps)[0],
            names=(self.names)[0]))

    def check_consistency(self):
        # check consistency of file list
        n_files = len(self.file_list)
        n_steps = len(self.steps)
        if n_files % n_steps != 0:
            raise RuntimeError('#files not divisible by #steps')
            # print('#files %i not divisible by #steps %i' % (n_files, n_steps))

        core_set = set()
        for name in self.names:
            for step in self.steps:
                new_list =  self.reduced_file_list(step, name)

                # get list of core numbers
                for filename in new_list:
                    z = re.match('.*' + name + '.+\.(.+)\.dat', filename)
                    if z:
                        core_set.add(int(z.groups()[0]))

                if not core_set:
                    raise RuntimeError('no files for ' + name + 'step ' + step)

                if max(core_set) != len(core_set) - 1 or min(core_set) != 0:
                    raise RuntimeError('file list inconsistent for ' + name
                            + ' step ' + step)

                for i in range(len(core_set) - 1):
                    if i not in core_set:
                        raise RuntimeError('no file for core no.' + str(i))

    def reduced_file_list(self, steps=None, names=None):
        if names is None:
            names = self.names
        else:
            if not isinstance(names, (list, tuple)):
                names = [names]

        if steps is None:
            steps = self.steps
        else:
            if not isinstance(steps, (list, tuple)):
                steps = [steps]

        new_list = []
        for step in steps:
            for name in names:
                for filename in self.file_list:
                    z = re.match('.*' + name + str(step) + '\.\d+\.dat', filename)

                    if z:
                        new_list.append(filename)

        new_list.sort()
        return new_list

    # copy function has to be defined individually in each subclass
    def copy(self, new_directory, steps=None, names=None):
        print(self.names, new_directory, steps, names)
        raise RuntimeError('function copy() not defined')

    def copy_first(self, new_directory, names=None):
        self.copy(new_directory, self.steps[0], names)

    def copy_last(self, new_directory, names=None):
        self.copy(new_directory, self.steps[-1], names)

    def info(self):
        print('Names:   ', self.names)
        print('#files = ', len(self.file_list))
        print('#steps = ', len(self.steps))


class Dump_folder(Dump_container):
    def __init__(self, path):
        self.path = path
        if not os.path.isdir(self.path):
            raise FileNotFoundError('file %s does not exist' % self.path)

        self.file_list = os.listdir(self.path)
        Dump_container.__init__(self, self.file_list)

    def copy(self, new_directory, steps=None, names=None):
        file_list = self.reduced_file_list(steps, names)
        if not file_list:
            print('No files to copy')
            return

        if not os.path.isdir(new_directory):
            print('creating directory ' + new_directory)
            os.mkdir(new_directory)


        for fl in file_list:
            copy2(self.path + '/' + fl, new_directory + '/' + fl)

    def concatenate(self, name, step, outfile):
        file_list = self.reduced_file_list(name, step)

        with open(outfile, 'wb') as wfd:
            for f in file_list:
                with open(f, 'rb') as fd:
                    #10MB per writing chunk to avoid reading big file into memory
                    copyfileobj(fd, wfd, 1024*1024*10)

    def compress(self, tarfile_name):
        try:
            tardir(self.path, tarfile_name)
        except:
            raise RuntimeError('tar command failed')


class Dump_tarfile(Dump_container):
    def __init__(self, path):
        self.path = path
        if not os.path.isfile(self.path):
            raise FileNotFoundError('file %s does not exist' % self.path)

        fl = []
        with tarfile.open(path, 'r') as tf:
            for tarinfo in tf:
                fl.append(tarinfo.name)

        self.file_list = sorted(fl)
        self.level = self.file_list[0].count('/')
        Dump_container.__init__(self, self.file_list)

    def reduced_member_list(self, tar, steps=None, names=None):
        if names is None:
            names = self.names
        else:
            if not isinstance(names, (list, tuple)):
                names = [names]

        if steps is None:
            steps = self.steps
        else:
            if not isinstance(steps, (list, tuple)):
                steps = [steps]

        new_list = []
        for step in steps:
            for name in names:
                for tarinfo in tar:
                    filename = tarinfo.name
                    if filename[0:2] == '..' or filename[0] == '/':
                        raise RuntimeError('Absolute or relative to parent '
                                'path detected\n%s' % filename)
                    z = re.match('.*' + name + str(step) + '\.\d+\.dat', filename)

                    if z:
                        new_list.append(tarinfo)

        return new_list

    def copy(self, new_directory, steps=None, names=None):
        if not os.path.isdir(new_directory):
            print('creating directory ' + new_directory)
            os.mkdir(new_directory)

        with tarfile.open(self.path, 'r') as tf:
            member_list = self.reduced_member_list(tf, steps, names)
            tf.extractall(new_directory, member_list)


def retrieve(directory, file_list):
    if not os.path.isdir(directory):
        raise FileNotFoundError('file %s does not exist' % directory)

    for fl in file_list:
        copy2(directory + '/' + fl, './' + fl)


def copy(directory, file_list, new_directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError('directory %s does not exist' % directory)

    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)

    for fl in file_list:
        copy2(directory + '/' + fl, new_directory + '/' + fl)


def cwd():
    ''' return current working directory with all links resolved '''
    wd = os.getcwd()
    return os.path.realpath(wd)


def scratch_path(path=cwd()):
    ''' return the path of the scratch directory that corresponds to path '''
    if not os.path.isdir(path):
        raise RuntimeError('Path {} does not exist'.format(path))

    USER = os.environ['USER']
    name = path

    # only take part relative to user directory
    name = re.sub('^.*/' + USER, '', name)
    # remove trailing slash(s)
    name = re.sub('/+$', '', name)
    # remove leading slash(s)
    name = re.sub('^/+', '', name)
    # replace remaining slashs by underscore
    name = re.sub('/', '_', name)

    # prepend scratch path
    if os.environ['TMPDIR']:
        prefix = os.environ['TMPDIR']
    else:
        prefix = '/ptmp/' + USER

    if not os.path.isdir(prefix):
        raise FileNotFoundError('directory ' + prefix + ' does not exist')

    name = prefix + '/' + name
    return name


def copy_to_scratch(src=cwd()):
    copy_tree(src, scratch_path())


def tardir(path, tar_name, mode='w'):
    with tarfile.open(tar_name, mode) as tar_handle:
        for root, _, files in os.walk(path):
            for f in files:
                tar_handle.add(os.path.join(root, f))


def archive_directory(directory, path, rm=False):
    print('archiving directory', directory)
    print('to path ', path)

    if not os.path.exists(path):
        print("creating " + path)
        os.mkdir(path)

    # archive rest of the data
    try:
        basename = os.path.abspath(directory)
        basename = os.path.basename(basename)

        temp_dir = tempfile.mkdtemp()
        temp_file = basename + '.tar'
        temp_path = temp_dir + '/' + temp_file

        archive_path = os.path.normpath(path + '/' + temp_file)
        if os.path.exists(archive_path):
            warnings.warn("WARNING: " + archive_path + ' already exists',
                    RuntimeWarning)
            return

        print('compressing directory to temporary file: ' + temp_path)
        returncode = subprocess.call(['tar', '-cf', temp_path, directory])
        if returncode != 0:
            raise ChildProcessError('tar command failed')

        print('moving archive')
        shutil.move(temp_path, path)

    finally:
        shutil.rmtree(temp_dir)

    if rm:
        print('removing ' + directory)
        shutil.rmtree(directory)


def archive_dump(directory, path, rm=False):
    print('archiving dump ' + directory)
    print('to path ', path)

    if not os.path.exists(path):
        print("creating " + path)
        os.mkdir(path)

    basename = os.path.abspath(directory)
    basename = os.path.dirname(basename)
    basename = os.path.basename(basename)

    try:
        temp_dir = tempfile.mkdtemp()
        temp_file = basename + '_dump.tar'
        temp_path = temp_dir + '/' + temp_file

        print('compressing dump to: ' + temp_path)
        returncode = subprocess.call(['tar', '-cf', temp_path, directory])
        if returncode != 0:
            raise ChildProcessError('tar command failed')

        archive_path = os.path.normpath(path + '/' + temp_file)
        if os.path.exists(archive_path):
            raise RuntimeError(archive_path + ' already exists')

        print('moving dump archive')
        shutil.move(temp_path, path)

        if rm and returncode == 0:
            print('removing ' + directory)
            shutil.rmtree(directory)

    finally:
        shutil.rmtree(temp_dir)
