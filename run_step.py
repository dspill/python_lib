import time
import os
import json
import espressopp as epp
from espressopp.tools.functions import setupSystem, customWritexyz, \
        customWritexyzStream, fileOutput, printInteractions, \
        arrangeChainsInSquareShape

def warmup(p):
    ''' Generate configuration and perform warmup '''
    seed               = 6543215 # seed for random
    temperature        = 1.0     # set temperature to None for NVE-simulations

    # set parameters for simulations
    num_chains = p['number_of_chains']
    monomers_per_chain = p['degree_of_polymerization']
    L = p['L']

    ######################################################################
    ### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
    ######################################################################

    nsteps      = p['num_steps_wu']
    isteps      = p['ints_per_step_wu']
    rc          = pow(2.0, 1.0/6.0)
    skin        = p['skin']
    timestep    = p['dt']
    box         = p['box']

    print(epp.Version().info())
    print('Setting up simulation ...')

    #logging.getLogger("SteepestDescent").setLevel(logging.INFO)

    system         = epp.System()
    system.rng     = epp.esutil.RNG()
    system.rng.seed(seed)
    system.bc      = epp.bc.OrthorhombicBC(system.rng, box)
    system.skin    = skin
    nodeGrid       = epp.tools.decomp.nodeGrid(epp.MPI.COMM_WORLD.size)
    cellGrid       = epp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
    system.storage = epp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    integrator     = epp.integrator.VelocityVerlet(system)
    integrator.dt  = timestep
    thermostat     = epp.integrator.LangevinThermostat(system)
    thermostat.gamma  = p['langevin_gamma']
    thermostat.temperature = p['temperature']
    integrator.addExtension(thermostat)

    steepest       = epp.integrator.MinimizeEnergy(system, gamma=0.001,
            ftol=0.1, max_displacement=0.001, variable_step_flag=False)

    # set the polymer properties
    bondlen            = 0.97

    props    = ['id', 'type', 'mass', 'pos', 'v']
    vel_zero = epp.Real3D(0.0, 0.0, 0.0)

    bondlist  = epp.FixedPairList(system.storage)
    #anglelist = epp.FixedTripleList(system.storage)
    pid       = 1
    bead_type = 0
    mass      = 1.0

    # add particles to the system and then decompose
    # do this in chunks of 1000 particles to speed it up
    chain = []
    for _ in range(num_chains):
        startpos = system.bc.getRandomPos()
        positions, bonds, _ = epp.tools.topology.polymerRW(pid,
                startpos, monomers_per_chain, bondlen, True)
        for k in range(monomers_per_chain):
            part = [pid + k, bead_type, mass, positions[k], vel_zero]
            chain.append(part)
        pid += monomers_per_chain
        #bead_type += 1
        system.storage.addParticles(chain, *props)
        system.storage.decompose()
        chain = []
        bondlist.addBonds(bonds)
        #anglelist.addTriples(angles)
    system.storage.addParticles(chain, *props)
    system.storage.decompose()

    num_particles = num_chains * monomers_per_chain
    density = num_particles * 1.0 / (L * L * L)

    # Lennard-Jones with Verlet list
    vl      = epp.VerletList(system, cutoff=rc)
    potLJ   = epp.interaction.LennardJones(epsilon=1.0, sigma=1.0,
            cutoff=rc, shift=0)
    interLJ = epp.interaction.VerletListLennardJones(vl)
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)

    # FENE bonds
    potFENE = epp.interaction.FENECapped(K=3000.0, r0=0.0, rMax=1.5,
                                                cutoff=8, r_cap=1.49999)
    interFENE = epp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
    system.addInteraction(interFENE, 'FENE')

    # Cosine with FixedTriple list
    #potCosine = epp.interaction.Cosine(K=1.5, theta0=3.1415926)
    #interCosine = epp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
    #system.addInteraction(interCosine)

    # print simulation parameters
    print('')
    print('number of particles = ', num_particles)
    print('length of system    = ', L)
    print('density             = ', density)
    print('rc                  = ', rc)
    print('dt                  = ', integrator.dt)
    print('skin                = ', system.skin)
    print('temperature         = ', temperature)
    print('nsteps              = ', nsteps)
    print('isteps              = ', isteps)
    print('NodeGrid            = ', system.storage.getNodeGrid())
    print('CellGrid            = ', system.storage.getCellGrid())
    print('')

    # epp.tools.decomp.tuneSkin(system, integrator)

    #filename = "initial_for_relax.res"
    #epp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)

    epp.tools.analyse.info(system, steepest)
    start_time = time.clock()
    for k in range(10):
        steepest.run(isteps)
        epp.tools.analyse.info(system, steepest)

    # exchange the FENE potential
    epp.System.removeInteractionByName(system, 'FENE')
    potFENE = epp.interaction.FENECapped(K=30.0, r0=0.0, rMax=1.5, cutoff=8, r_cap=1.499999)
    interFENE = epp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
    system.addInteraction(interFENE)

    for k in range(20):
        steepest.run(isteps)
        epp.tools.analyse.info(system, steepest)
    end_time = time.clock()

    epp.tools.analyse.info(system, integrator)
    for k in range(2):
        integrator.run(isteps)
        epp.tools.analyse.info(system, integrator)
    end_time = time.clock()
    epp.tools.analyse.info(system, integrator)
    epp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

    customWritexyz("final_configuration_warmup.xyz", system, velocities=True,
                   append=False)


def warmup2d(p):
    ''' Generate configuration, put it in the z=0 plane and perform warmup '''
    timeout = time.time() + 0.97 * p['timelimit'] - 90
    simstep = 'warmup2d'

    system, integrator, _ = setupSystem(p, xyzfilename=None, phi=0.,
                                         with_lb=False)
    printInteractions(system)

    arrangeChainsInSquareShape(system, p['number_of_chains'],
            p['degree_of_polymerization'], p['box'][0])

    # display observables
    epp.tools.analyse.info(system, integrator)
    # output observables
    outname = "output_"+simstep+".dat"
    fileOutput(system, integrator, outname)
    # output trajectory
    trajname = "traj_"+simstep+".xyz"
    filestream = open(trajname, 'w')
    customWritexyzStream(filestream, system)

    duration = 0.
    while time.time() < (timeout - duration):
        duration = time.time()

        integrator.run(p['ints_per_step'])
        # display observables
        epp.tools.analyse.info(system, integrator)
        # output observables
        fileOutput(system, integrator, outname)
        # output trajectory
        customWritexyzStream(filestream, system)

        duration = time.time() - duration

    filestream.close()

    customWritexyz('final_configuration_' + simstep + '.xyz', system,
            velocities=True,
            append=False)


def run(simstep, p, xyzfilename):
    ''' Run simulation with parameters p '''
    if simstep == 'quench':
        fac = 0.9
    else:
        fac = 0.97

    timeout = time.time() + fac * p['timelimit'] - 90

    print("Running %s" % simstep)
    if simstep == 'relax':
        phi = 0.
        with_lb = False
    elif simstep == 'add_lb':
        phi = 0.
        with_lb = True
    elif simstep == 'quench':
        phi = p['phi']
        with_lb = True
    else:
        raise RuntimeError("simstep " + simstep + " not supported")

    if not os.path.isfile(xyzfilename):
        raise RuntimeError("File %s does not exist" % xyzfilename)

    system, integrator, lb = setupSystem(p, xyzfilename=xyzfilename, phi=phi,
                                         with_lb=with_lb)
    printInteractions(system)

    # display observables
    epp.tools.analyse.info(system, integrator)
    # output observables
    outname = "output_"+simstep+".dat"
    fileOutput(system, integrator, outname)
    # output trajectory
    trajname = "traj_"+simstep+".xyz"
    filestream = open(trajname, 'w')
    customWritexyzStream(filestream, system)

    duration = 0.
    while time.time() < (timeout - duration):
        duration = time.time()

        integrator.run(p['ints_per_step'])
        # display observables
        epp.tools.analyse.info(system, integrator)
        # output observables
        fileOutput(system, integrator, outname)
        # output trajectory
        customWritexyzStream(filestream, system)

        if lb:
            # output LB configuration
            if simstep != 'add_lb':
                lb.keepLBDump()     # flag to keep previously saved LB state

            lb.saveLBConf()     # saves current state of the LB fluid

        duration = time.time() - duration

    filestream.close()

    customWritexyz('final_configuration_' + simstep + '.xyz', system,
                   velocities=True, append=False)


def read_parameters(filename):
    p = json.load(open(filename))
    p['box'] = tuple(p['box'])
    return p
