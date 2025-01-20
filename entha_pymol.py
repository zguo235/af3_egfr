# %%
from simtk.openmm import app
from simtk.openmm import *
from simtk.unit import *
import pymol
from pymol import cmd

# %%
# Initialize PyMOL
pymol.finish_launching(['pymol', '-cq'])

# Load the CIF file
cmd.load("/home/guo1amed/Projects/AF3/af_output/protein_complex_1jrh/seed-1_sample-0/model.cif", "model")

# Add hydrogens
cmd.h_add("model")

# Save as PDB file
cmd.save("input.pdb", "model")

# %%
# Load the PDB file
pdb = app.PDBFile('input.pdb')

# Create a force field
forcefield = app.ForceField('amber99sb.xml', 'tip3p.xml')

# Create a system
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=1*nanometer, constraints=app.HBonds)

# Create an integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Create a simulation
simulation = app.Simulation(pdb.topology, system, integrator)

# Set the initial positions
simulation.context.setPositions(pdb.positions)

# Minimize the energy
simulation.minimizeEnergy()

# Get the state of the system
state = simulation.context.getState(getEnergy=True)

# Calculate the enthalpy
potential_energy = state.getPotentialEnergy()
kinetic_energy = state.getKineticEnergy()
enthalpy = potential_energy + kinetic_energy

print(f'Potential Energy: {potential_energy}')
print(f'Kinetic Energy: {kinetic_energy}')
print(f'Enthalpy: {enthalpy}')