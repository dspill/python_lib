import os
import json
import datetime
import argparse
from shutil import rmtree
from functions import hours_minutes_seconds, YesNo
from file_operations import scratch_path, cwd, Dump_folder

valid_simsteps = ['warmup', 'warmup2d', 'relax', 'add_lb', 'quench', 'cleanup']

class Submit_script:
    # def __init__(self, **parameters):
    def __init__(self, parameters, force=False):
        self.parameters        = parameters

        self.simstep           = parameters['simstep']
        self.partition         = parameters['partition']
        self.n_nodes           = parameters['n_nodes']
        self.n_tasks_per_node  = parameters['n_tasks_per_node']
        self.n_tasks_per_core  = parameters['n_tasks_per_core']
        self.timelimit         = parameters['timelimit']
        self.email             = parameters['email']
        self.mail_type         = parameters['mail_type']
        self.esprc_path        = parameters['esprc_path']
        self.infile            = parameters['infile']
        self.machine           = parameters['machine']
        self.lattice_boltzmann = parameters['lattice_boltzmann']

        self.name = self.simstep
        if 'name' in parameters:
            if parameters['name'] is not None:
                self.name = parameters['name']

        self.force             = force
        self.para_filename     = 'parameters_' + self.simstep + '.txt'
        self.n_cores           = self.n_nodes * self.n_tasks_per_node
        self.n_tasks           = self.n_cores * self.n_tasks_per_core

        print('  Running ' + str(self.n_tasks) + ' tasks on '
                + str(self.n_cores) + ' cores')
        print('  Timelimit: '
                + str(datetime.timedelta(seconds=self.timelimit)))

        if not force:
            self.check_consistency()

    def check_infile(self):
        # check parameters of infile
        f = open(self.infile)

        line = f.readline()
        number_of_particles = int(line.split()[0])
        line = f.readline().split()
        if len(line) == 3:
            Lx = float(line[0])
            Ly = float(line[1])
            Lz = float(line[2])
        else:
            Lx = float(line[0])
            Ly = float(line[4])
            Lz = float(line[8])
        box = (Lx, Ly, Lz)

        if number_of_particles != self.parameters['number_of_particles']:
            raise RuntimeError('num particles from input file does not match')
        if box != self.parameters['box']:
            raise RuntimeError('box from input file does not match')

    def check_consistency(self):
        variables = {
            'simstep'          : self.simstep,
            'partition'        : self.partition,
            'n_nodes'          : self.n_nodes,
            'n_tasks_per_node' : self.n_tasks_per_node,
            'n_tasks_per_core' : self.n_tasks_per_core,
            'timelimit'        : self.timelimit,
            'email'            : self.email,
            'mail_type'        : self.mail_type,
            'esprc_path'       : self.esprc_path,
            'machine'          : self.machine,
            'lattice_boltzmann': self.lattice_boltzmann,
            'name'             : self.name,
            }
        for key, var in variables.items():
            if var is None:
                raise RuntimeError('Variable ' + key + ' not set')

        if self.simstep not in valid_simsteps:
            raise RuntimeError('simstep ' + self.simstep + ' not supported')

        if self.simstep not in ['warmup', 'warmup2d', 'cleanup']:
            infile = self.infile
            if not infile:
                raise RuntimeError('You have to give an input file')

            infile = os.path.expanduser(infile)
            infile = os.path.normpath(infile)
            if not os.path.isfile(infile):
                raise FileNotFoundError(infile + ' does not exist')

            self.check_infile()

        self.esprc_path = os.path.expanduser(self.esprc_path)
        self.esprc_path = os.path.normpath(self.esprc_path)
        if not os.path.isfile(self.esprc_path):
            raise FileNotFoundError(self.esprc_path + ' does not exist')

    def write_slurm(self, outfile):
        ''' write slurm directives '''
        outfile.write('#!/bin/bash\n')
        outfile.write('# Standard output and error:\n')
        outfile.write('#SBATCH -o ./tjob_' + self.simstep + '.out.%j\n') # log stdout
        outfile.write('#SBATCH -e ./tjob_' + self.simstep + '.err.%j\n') # log stderr
        outfile.write('\n')

        outfile.write('#SBATCH -D ./\n') # working directory
        outfile.write('#SBATCH -J ' + self.name + '\n') # job name
        outfile.write('\n')

        outfile.write('#SBATCH --partition=' + self.partition + '\n')
        if self.partition not in ['express', 'small']:
            outfile.write('#SBATCH --nodes=' + str(self.n_nodes)+'\n')
        outfile.write('#SBATCH --ntasks-per-node='+str(self.n_tasks_per_node)+'\n')
        outfile.write('#SBATCH --ntasks-per-core='+str(self.n_tasks_per_core)+'\n')
        timelimit_datetime = datetime.timedelta(seconds=self.timelimit)
        outfile.write('#SBATCH --time='
                + hours_minutes_seconds(timelimit_datetime))
        outfile.write('\n')

        if self.email:
            outfile.write('#SBATCH --mail-type=' + self.mail_type + '\n')
        outfile.write('#SBATCH --mail-user=' + self.email + '\n')
        outfile.write('\n')

    def write_modules(self, outfile):
        outfile.write('module purge\n')

        if self.machine == 'draco':
            outfile.write('module load intel/18.0.5\n')
            outfile.write('module load impi/2018.4\n')
            outfile.write('module load anaconda/3/5.1\n')

        elif self.machine == 'cobra':
            outfile.write('module load intel\n')
            outfile.write('module load impi/2018.4\n')
            outfile.write('module load mkl\n')
            outfile.write('module load anaconda/3/5.1\n')

        else:
            raise RuntimeError(self.machine + 'is not a known machine')

        outfile.write('\n')

    def write_run_code(self, outfile):
        outfile.write('srun python2.7 << END\n')
        outfile.write('import run_step\n')
        outfile.write('p = run_step.read_parameters(''\''
                + self.para_filename + '\')\n')
        if self.simstep == 'warmup':
            outfile.write('run_step.warmup(p)\n')
        elif self.simstep == 'warmup2d':
            outfile.write('run_step.warmup2d(p)\n')
        else:
            outfile.write('run_step.run(\'' + self.simstep
                    + '\', p, \'' + os.path.basename(self.infile) + '\')\n')

        outfile.write('END\n')


    def write_submit_script(self, scriptname):
        with open(scriptname, 'w') as outfile:
            self.write_slurm(outfile)

            outfile.write('echo "Job $SLURM_JOB_ID"\n')

            outfile.write('\n')
            self.write_modules(outfile)

            # setup python environment
            outfile.write('PYTHONPATH=""\n')
            if not os.path.isfile(os.environ['HOME'] + '/python/run_step.py'):
                raise RuntimeError('File ' + os.environ['HOME']
                        + '/python/run_step.py does not exist')

            outfile.write('PYTHONPATH="$HOME/python"\n')

            outfile.write('source ' + self.esprc_path)
            outfile.write('\n')

            # set some variables
            scratch = scratch_path()
            if os.path.isdir(scratch) and self.simstep != 'cleanup':
                if YesNo('Scratch directory ' + scratch
                        + ' exists. Delete it? '):
                    print('deleting...')
                    rmtree(scratch)
                    print('...done')


            wd = cwd()
            outfile.write('wd=' + wd + '\n')
            outfile.write('scratch=' + scratch + '\n')
            outfile.write('\n')

            # clean scratch_path
            outfile.write('echo \"clearing scratch folder\"\n')
            outfile.write('rm -rf $scratch\n')
            outfile.write('mkdir -pv $scratch\n')
            # copy scripts to scratch
            outfile.write('rsync -rvat *.{sh,txt} $scratch 2>&1\n')

            # move initial configurations to scratch
            if self.infile:
                outfile.write('rsync -rvat ' + self.infile
                        + ' $scratch 2>&1\n')
                infile_base = os.path.dirname(self.infile)


            if self.lattice_boltzmann:
                # look in current directory for dump folder
                if os.path.isdir('./dump'):
                    initial_dump = './dump'
                # look in directory of input file for dump folder
                elif self.infile and os.path.isdir(infile_base + '/dump'):
                    initial_dump = os.path.normpath(infile_base + '/dump')
                # for the quench we require initial dump
                elif not self.force and self.simstep == 'quench':
                    raise RuntimeError('No initial LB configuration found')
                else:
                    initial_dump = None

                if initial_dump is None:
                    print('starting with new LB configuration')
                else:
                    print('using most recent configuration in %s as '
                    'initial LB configuration' % initial_dump)
                    df = Dump_folder(initial_dump)
                    if df.n_cores() != self.n_cores:
                        raise RuntimeError("number of cores does not match")
                    # copy most recent LB configuration to scratch
                    copy_script='~/scripts/copy_last_LB_configuration.py'
                    if not os.path.isfile(os.path.expanduser(copy_script)):
                        raise FileNotFoundError(copy_script + ' does not exist')
                    outfile.write(copy_script + ' ' + initial_dump + ' $scratch/dump\n')

            outfile.write('\n')
            outfile.write('cd $scratch\n')

            # Run the program:
            self.write_run_code(outfile)
            outfile.write('\n')
            outfile.write('result=$?\n')
            outfile.write('\n')

            outfile.write('if [ $result -eq 0 ]; then\n')
            outfile.write('    echo \'job successful\'\n')
            # retrieve files
            outfile.write('    rsync -rvat '
                    '*.{dat,txt} final_configuration_*.xyz $wd 2>&1\n')
            if self.lattice_boltzmann:
                copy_script='~/scripts/copy_last_LB_configuration.py'
                if not os.path.isfile(os.path.expanduser(copy_script)):
                    raise FileNotFoundError(copy_script + ' does not exist')
                outfile.write('    ' + copy_script + ' $scratch/dump $wd/dump\n')
            outfile.write('    echo DONE\n')
            outfile.write('    exit $result\n')

            outfile.write('else\n')
            outfile.write('    echo \'job FAILED\'\n')
            outfile.write('    echo DONE\n')
            outfile.write('    exit $result\n')
            outfile.write('fi\n')

            outfile.write('echo DONE\n')

    def write_parameters(self):
        json.dump(self.parameters, open(self.para_filename,'w'))

    def write(self, submit_scriptname):
        self.write_submit_script(submit_scriptname)
        self.write_parameters()


    def write_cleanup_script(self, scriptname, backupdir):
        with open(scriptname, 'w') as outfile:
            self.write_slurm(outfile)
            outfile.write('echo "Job $SLURM_JOB_ID"\n')
            outfile.write('\n')

            # set some variables
            scratch = scratch_path()

            # Run the program:
            outfile.write('python3 << END\n')
            outfile.write('import file_operations as fo\n')
            outfile.write('fo.archive_dump(\'' + scratch + '/dump\', ' 
                    + '\'' + backupdir + '\'' + ', rm=True)\n')
            outfile.write('END\n')

            outfile.write('\n')
            outfile.write('result=$?\n')
            outfile.write('\n')

            outfile.write('if [ $result -eq 0 ]; then\n')
            outfile.write('    echo \'job successful\'\n')
            outfile.write('    echo DONE\n')
            outfile.write('    exit $result\n')

            outfile.write('else\n')
            outfile.write('    echo \'job FAILED\'\n')
            outfile.write('    echo DONE\n')
            outfile.write('    exit $result\n')
            outfile.write('fi\n')

            outfile.write('echo DONE\n')

def generate(p):
    ''' generate a submit script that can be submitted via the command
    sbatch
    '''
    arguments = parseArguments()
    p['simstep']   = arguments['simstep']
    p['infile']    = arguments['infile']
    p['name']      = arguments['name']

    scriptname = 'submit_' + p['simstep'] + '.sh'

    print('Generating ' + scriptname)
    # warmup on express partition
    if p['simstep'] in ['warmup', 'warmup2d']:
        p['partition']         = 'express'
        p['n_nodes']           = 1
        p['n_tasks_per_node']  = 32
        p['n_tasks_per_core']  = 1
        p['timelimit']         = 30*60
        p['lattice_boltzmann'] = False
    elif p['simstep'] == 'cleanup':
        p['partition']         = 'small'
        p['n_nodes']           = 1
        p['n_tasks_per_node']  = 1
        p['n_tasks_per_core']  = 1
        p['timelimit']         = datetime.timedelta(hours=12).total_seconds()
        p['lattice_boltzmann'] = False
    elif p['simstep'] == 'relax':
        p['lattice_boltzmann'] = False
    elif p['simstep'] == 'add_lb':
        p['lattice_boltzmann'] = True
    elif p['simstep'] == 'quench':
        p['lattice_boltzmann'] = True
    else:
        raise RuntimeError('simstep ' + p['simstep'] + ' not valid')

    # ================== TESTING ======================
    if arguments['test']:
        p['partition']         = 'express'
        p['n_nodes']           = 1
        p['n_tasks_per_node']  = 16
        p['n_tasks_per_core']  = 1
        p['timelimit']         = 5*60
        p['ints_per_step']     = 100
    # ================== TESTING ======================


    submit_script = Submit_script(p, force=arguments['force'])
    if arguments['simstep'] == 'cleanup':
        submit_script.write_cleanup_script(scriptname, arguments['backupdir'])
    else:
        # write parameters to file
        submit_script.write_parameters()
        # write actual submit script
        submit_script.write(scriptname)


def parseArguments():
    parser = argparse.ArgumentParser(
        description='''Generate slurm submit script''')
    parser.add_argument('--simstep', '-s', type=str,
            choices=valid_simsteps, required=True)
    parser.add_argument('--infile', '-i', type=str, help='input file')
    parser.add_argument('--test', '-t', action='store_true', help='Short test run')
    parser.add_argument('--force', '-f', action='store_true',
            help="Force generation and do not sanity check")
    parser.add_argument('--name', '-n', type=str,
            help='name: slurm job name', required=False, default=None)
    parser.add_argument('--backupdir', '-b', type=str, help='backup directory',
            required=False, default='/r/d/dspiller/')

    arguments = vars(parser.parse_args())
    if arguments['simstep'] == 'cleanup' and not arguments['backupdir']:
        raise RuntimeError('you have to specify --backupdir')
    return arguments
