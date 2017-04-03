#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for a Polymer Melt System including           #
#  runtime details
#                                                                         #
###########################################################################

import time
import espressopp

nsteps      = 10
isteps      = 100
rc          = pow(2.0, 1.0/6.0)
skin        = 0.4
timestep    = 0.001
density     = 0.27
num_particles = 6000
L = pow(num_particles/density, 1.0/3.0)
box = (L, L, L)
# FENE
r0 = 0.97
K = 30.0
# Angles
kappa = 2.5
theta0 = 3.1415
capradius          = 0.6
equil_nloops       = 100
equil_isteps       = 100
warmup_cutoff      = pow(2.0, 1.0/6.0)
warmup_nloops      = 100
warmup_isteps      = 200
total_warmup_steps = warmup_nloops * warmup_isteps
epsilon_start      = 0.1
epsilon_end        = 1.0
epsilon_delta      = (epsilon_end - epsilon_start) / warmup_nloops
epsilon            = 1.0
sigma              = 1.0

# set temperature to None for NVE-simulations
temperature = 1.0

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print espressopp.Version().info()
print 'Setting up simulation ...'
system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

# Create system
angles = [(i+1, i+2, i+3) for i in range(0, num_particles, 3)]
bonds = [(b[0], b[1]) for b in angles] + [(b[1], b[2]) for b in angles]

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(0, num_particles, 3):
  pos = system.bc.getRandomPos()
  new_particles.append([i+1, 0, 1.0, pos])
  new_particles.append([i+2, 0, 1.0, pos+espressopp.Real3D(0, 0, r0)])
  new_particles.append([i+3, 0, 1.0, pos+espressopp.Real3D(0, 0, 2*r0)])
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)
system.storage.decompose()


# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espressopp.interaction.Harmonic(K=30.0, r0=r0)
interFENE = espressopp.interaction.FixedPairListHarmonic(system, fpl, potFENE)
system.addInteraction(interFENE, 'fene')

# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espressopp.interaction.Cosine(K=kappa, theta0=theta0)
interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
system.addInteraction(interCosine, 'cosine')

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

########################################################################
# 5. setting up interaction potential for the warmup                   #
########################################################################

# create a verlet list that uses a cutoff radius = warmup_cutoff
# the verlet radius is automatically increased by system.skin (see system setup)
verletlist  = espressopp.VerletList(system, rc)
# create a force capped Lennard-Jones potential
# the potential is automatically shifted so that U(r=cutoff) = 0.0
LJpot       = espressopp.interaction.LennardJonesCapped(epsilon=epsilon_start, sigma=sigma, cutoff=rc, caprad=capradius, shift='auto')
# create a force capped Lennard-Jones interaction that uses a verlet list 
interaction = espressopp.interaction.VerletListLennardJonesCapped(verletlist)
# tell the interaction to use the above defined force capped Lennard-Jones potential
# between 2 particles of type 0 
interaction.setPotential(type1=0, type2=0, potential=LJpot)

########################################################################
# 6. running the warmup loop
########################################################################

# make the force capping interaction known to the system
system.addInteraction(interaction, 'ljcapp')
print "starting warm-up ...", warmup_nloops*warmup_isteps
# print some status information (time, measured temperature, pressure,
# pressure tensor (xy only), kinetic energy, potential energy, total energy, boxsize)
espressopp.tools.analyse.info(system, integrator)
for step in range(warmup_nloops):
  # perform warmup_isteps integraton steps
  integrator.run(warmup_isteps)
  # decrease force capping radius in the potential
  LJpot.epsilon += epsilon_delta
  # update the type0-type0 interaction to use the new values of LJpot
  interaction.setPotential(type1=0, type2=0, potential=LJpot)
  # print status info
  espressopp.tools.analyse.info(system, integrator)  
print "warmup finished"

system.removeInteractionByName('ljcapp')

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift='auto')
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# espressopp.tools.decomp.tuneSkin(system, integrator)

dump_gro = espressopp.io.DumpGRO(system, integrator, 'dump.gro', unfolded=True, append=True)
extanal = espressopp.integrator.ExtAnalyze(dump_gro, 10)
integrator.addExtension(extanal)

espressopp.tools.analyse.info(system, integrator)
start_time = time.clock()
for k in range(nsteps):
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)
end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

dump_final = espressopp.io.DumpGRO(system, integrator, 'final.gro', unfolded=True)
dump_final.dump()
