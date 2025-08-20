from argparse import ArgumentParser
import math
import numpy as np
import matplotlib.pyplot as plt

import openmc
import sys
import os
import csv
from array import array

# Get common input parameters shared by other physics
script_dir = os.path.dirname(__file__)
sys.path.append(script_dir)
import common_input as specs

model = openmc.model.Model()
# -------------- Unit Conversions: OpenMC requires cm -----------
m = 100.0

# geometry
reactor_height = 0.025*m     #April specified that number
reactor_bottom = -reactor_height/2.0
reactor_top = reactor_bottom + reactor_height

# geometry parameters
fuel_channel_diam = specs.compact_diameter * m


# superimposed search lattice
#triso_lattice_shape = (4, 4, 1)
triso_lattice_shape = (1, 1, 1)

all_materials = []

enrichment = 0.155        # U-235 enrichment (weight percent)
enrichment_234 = 2e-3     # U-234 enrichment (weight percent)
kernel_density = 10820    # fissile kernel density (kg/m3)
buffer_density = 1050     # buffer density (kg/m3)
PyC_density = 1900        # PyC density (kg/m3)
SiC_density = 3203        # SiC density (kg/m3)
matrix_density = 1700     # graphite matrix density (kg/m3)
helium_density = 0.18     # graphite matrix density (kg/m3)

# ----- uranium oxycarbide fuel ----- #
m_fuel = openmc.Material(name='fuel')
mass_234 = openmc.data.atomic_mass('U234')
mass_235 = openmc.data.atomic_mass('U235')
mass_238 = openmc.data.atomic_mass('U238')

n_234 = enrichment_234 / mass_234
n_235 = enrichment / mass_235
n_238 = (1.0 - enrichment - enrichment_234) / mass_238
total_n = n_234 + n_235 + n_238

m_fuel.add_nuclide('U234', n_234 / total_n)
m_fuel.add_nuclide('U235', n_235 / total_n)
m_fuel.add_nuclide('U238', n_238 / total_n)
m_fuel.add_element('C'   , 1.50)
m_fuel.add_element('O'   , 0.50)
m_fuel.set_density('kg/m3', kernel_density)
all_materials.append(m_fuel)


# ----- graphite buffer ----- #
m_graphite_c_buffer = openmc.Material(name='buffer')
m_graphite_c_buffer.add_element('C', 1.0)
m_graphite_c_buffer.add_s_alpha_beta('c_Graphite')
m_graphite_c_buffer.set_density('kg/m3', buffer_density)
all_materials.append(m_graphite_c_buffer)


# ----- pyrolitic carbon ----- #
m_graphite_pyc = openmc.Material(name='pyc')
m_graphite_pyc.add_element('C', 1.0)
m_graphite_pyc.add_s_alpha_beta('c_Graphite')
m_graphite_pyc.set_density('kg/m3', PyC_density)
all_materials.append(m_graphite_pyc)

# ----- silicon carbide ----- #
m_sic = openmc.Material(name='sic')
m_sic.add_element('C' , 1.0)
m_sic.add_element('Si', 1.0)
m_sic.set_density('kg/m3', SiC_density)
all_materials.append(m_sic)

# ----- matrix graphite ----- #
m_graphite_matrix = openmc.Material(name='graphite moderator')
m_graphite_matrix.add_element('C', 1.0)
m_graphite_matrix.add_s_alpha_beta('c_Graphite')
m_graphite_matrix.set_density('kg/m3', matrix_density)
all_materials.append(m_graphite_matrix)

# Coolant
m_coolant = openmc.Material(name='Helium coolant')
m_coolant.add_element('He', 1.0, 'ao')
m_coolant.temperature = 760.0
m_coolant.set_density('kg/m3', helium_density)
all_materials.append(m_coolant)  


materials = openmc.Materials(all_materials)
materials.export_to_xml()

# TRISO particle
radius_pyc_outer  = specs.oPyC_radius * m
radius_graphite   = specs.graphite_radius * m

s_fuel             = openmc.Sphere(r=specs.kernel_radius*m)
s_c_buffer         = openmc.Sphere(r=specs.buffer_radius*m)
s_pyc_inner        = openmc.Sphere(r=specs.iPyC_radius*m)
s_sic              = openmc.Sphere(r=specs.SiC_radius*m)
s_pyc_outer        = openmc.Sphere(r=radius_pyc_outer)
s_graphite         = openmc.Sphere(r=radius_graphite)
c_triso_fuel       = openmc.Cell(name='c_triso_fuel'     , fill=m_fuel,              region=-s_fuel)
c_triso_c_buffer   = openmc.Cell(name='c_triso_c_buffer' , fill=m_graphite_c_buffer, region=+s_fuel      & -s_c_buffer)
c_triso_pyc_inner  = openmc.Cell(name='c_triso_pyc_inner', fill=m_graphite_pyc,      region=+s_c_buffer  & -s_pyc_inner)
c_triso_sic        = openmc.Cell(name='c_triso_sic'      , fill=m_sic,               region=+s_pyc_inner & -s_sic)
c_triso_pyc_outer  = openmc.Cell(name='c_triso_pyc_outer', fill=m_graphite_pyc,      region=+s_sic       & -s_pyc_outer)
c_triso_matrix     = openmc.Cell(name='c_triso_matrix'   , fill=m_graphite_matrix,   region=+s_pyc_outer)
u_triso            = openmc.Universe(cells=[c_triso_fuel, c_triso_c_buffer, c_triso_pyc_inner, c_triso_sic, c_triso_pyc_outer, c_triso_matrix])

# Channel Fuel surface
fuel_cyl = openmc.ZCylinder(r = 0.5 * 1.27) # Radius of the fuel cylinder (1.27 diameter)
coolant_cyl = openmc.ZCylinder(r = 0.5 * 1.588)  # Coolant channel radius (1.588 diameter)


# create a TRISO lattice 
min_z = openmc.ZPlane(z0 = reactor_bottom,  boundary_type = 'reflective')
max_z = openmc.ZPlane(z0 = reactor_top,  boundary_type = 'reflective')

# region in which TRISOs are generated
r_triso = -fuel_cyl & +min_z & -max_z

rand_spheres = openmc.model.pack_spheres(radius=radius_pyc_outer , region=r_triso, pf=0.4, seed=1)
# Cartesian coordinates of sphere centers

# Export the coordinates of the sphere centers to csv file
with open('output.csv', 'w', newline='') as f_output:
    csv_output = csv.writer(f_output)
    csv_output.writerows(rand_spheres * 1e4)
    
random_trisos = [openmc.model.TRISO(radius_pyc_outer  , u_triso, i) for i in rand_spheres]

llc, urc = r_triso.bounding_box
pitch = (urc - llc) / triso_lattice_shape

# insert TRISOs into a lattice to accelerate point location queries
triso_lattice = openmc.model.create_triso_lattice(random_trisos, llc, pitch, triso_lattice_shape, m_graphite_matrix)

# Create a universe to contain the TRISO lattice
triso_universe = openmc.Universe(name='TRISO Universe')
triso_universe.add_cell(openmc.Cell(fill=triso_lattice))

# Number of cylinders
n_c = 10

# Total radius of the original cylinder
total_radius = 0.5 * 1.27  # 1.27 diameter divided by 2 for radius

# Create the smaller cylindrical surfaces
cylindrical_surfaces = []
for i in range(n_c):
    #r = (i + 1) * total_radius / n_c
    r = math.sqrt((i + 1) / n_c) * total_radius 
    cylindrical_surfaces.append(openmc.ZCylinder(r=r))

# Create the cells for each smaller cylinder
fuel_cells = []
for i in range(n_c):
    if i == 0:
        region = -cylindrical_surfaces[i]
    else:
        region = +cylindrical_surfaces[i - 1] & -cylindrical_surfaces[i]
    fuel_cells.append(openmc.Cell(region=region, fill=triso_universe, name=f'fuel_cell_{i + 1}'))


print(fuel_cells)
# Create the outer cell (graphite matrix surrounding the fuel)
outer_region = +cylindrical_surfaces[-1]
outer_cell = openmc.Cell(region=outer_region, fill=m_graphite_matrix)

# Create a universe for the fuel
fuel_universe = openmc.Universe(name='fuel_universe')
fuel_universe.add_cells(fuel_cells)
fuel_universe.add_cell(outer_cell)
print(fuel_universe)

coolant_cell = openmc.Cell(region=-coolant_cyl, fill=m_coolant)  # coolant cylinder 
coolant_ch_matrix_cell = openmc.Cell(region=+coolant_cyl, fill=m_graphite_matrix)
coolant_universe = openmc.Universe(cells=[coolant_cell, coolant_ch_matrix_cell])

all_graphite_cell = openmc.Cell(fill=m_graphite_matrix)
outer_universe = openmc.Universe(cells=[all_graphite_cell])

# None universe for the hexagonal lattice
none_cell = openmc.Cell(region=-coolant_cyl, fill=None)  # coolant cylinder 
none_ch_matrix_cell = openmc.Cell(region=+coolant_cyl, fill=None)
none_universe = openmc.Universe(cells=[none_cell, none_ch_matrix_cell])

lattice = openmc.HexLattice()
lattice.center = (0., 0.)
lattice.pitch = (1.8799,)  # Distance between the center of ring 1 to ring 2, calculated using the distance between two coolant channels, it's 3.256/sqrt(3)

lattice.outer = outer_universe
lattice.orientation = 'x'  # To match the planes of the heat conduction mesh
print(lattice.show_indices(num_rings=2))

#outer_ring = [fuel_universe, coolant_universe, fuel_universe, coolant_universe, fuel_universe, coolant_universe]  # Adds up to 6
outer_ring = [none_universe, coolant_universe, none_universe, coolant_universe, none_universe, coolant_universe]  # Adds up to 6
inner_ring = [fuel_universe]

lattice.universes = [outer_ring, inner_ring]
print(lattice)

# 3 planes, using the points (0.9385, 1.62553, 1.25), (0.9385, 1.62553, -1.25), (0.9385, -1.62553, 1.25), (0.9385, -1.62553, -1.25), (1.8799, 0, 1.25), and (1.8799, 0, -1.25)
plane1 = openmc.Plane(a=3.251, b=0.0, c=0.0, d=3.0510635, boundary_type='reflective')
plane2 = openmc.Plane(a=1.6255, b=-2.8184, c=0.0, d=-3.05577745, boundary_type='reflective')
plane3 = openmc.Plane(a=-1.6255, b=-2.8184, c=0.0, d=3.05577745, boundary_type='reflective')

main_cell = openmc.Cell(fill=lattice, region=+min_z & -max_z & -plane1 & +plane2 & -plane3)
geometry = openmc.Geometry([main_cell])

geometry.export_to_xml()

# Bounding box for calculating the fuel volume
lower_left_volume = (-2.0, -2.0, -1.25)
upper_right_volume = (2.0, 2.0, 1.25)
vol_calc = openmc.VolumeCalculation([m_fuel], 1000000, lower_left_volume, upper_right_volume)

### Settings ######
settings = openmc.Settings()

n_inactive = 500
n_active = 1500
settings.particles = 10000
settings.inactive = n_inactive
settings.batches = settings.inactive + n_active
settings.temperature = {'default': 760.0,
                        'method': 'interpolation',
                        'multipole': True,
                        'range': (0.0, 3000.0)}
lower_left = (-1.0, -1.0, reactor_bottom)
upper_right = (1.0, 1.0, reactor_top)
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)
settings.volume_calculations = [vol_calc]

# Model settings
model.settings = settings
settings.export_to_xml()

cell_pitch = 3.256

plot1 = openmc.Plot()
plot1.filename = 'plot1'
plot1.width = (3 * cell_pitch, 3 * cell_pitch)
plot1.basis = 'xz'
plot1.origin = (0.0, 0.0, 0.0)
plot1.pixels = (int(800 * 3.256), int(800 * 3.256))
plot1.color_by = 'cell'
plot1.level = 1

plot5 = openmc.Plot()
plot5.filename = 'plot5'
plot5.width = (3 * cell_pitch, 3 * cell_pitch)
plot5.basis = 'xz'
plot5.origin = (0.0, 0.0, 0.0)
plot5.pixels = (int(800 * 3.256), int(800 * 3.256))
plot5.color_by = 'cell'


plot2 = openmc.Plot()
plot2.filename = 'plot2'
plot2.width = (3 * cell_pitch, 3 * cell_pitch)
plot2.basis = 'xy'
plot2.origin = (0.0, 0.0, 0.0)
plot2.pixels = (int(800 * 3.256), int(800 * 3.256))
plot2.color_by = 'material'


plot3 = openmc.Plot()
plot3.filename = 'plot3'
plot3.width = plot2.width
plot3.basis = plot2.basis
plot3.origin = plot2.origin
plot3.pixels = plot2.pixels
plot3.color_by = 'cell'
plot3.level = 1

plot6 = openmc.Plot()
plot6.filename = 'plot6'
plot6.width = plot2.width
plot6.basis = plot2.basis
plot6.origin = plot2.origin
plot6.pixels = plot2.pixels
plot6.color_by = 'cell'


plot4 = openmc.Plot()
plot4.filename = 'plot4'
plot4.width = (3 * cell_pitch, 3 * cell_pitch)
plot4.basis = 'yz'
plot4.origin = (0.0, 0.0, 0.0)
plot4.pixels = (int(800 * 3.256), int(800 * 3.256))
plot4.color_by = 'material'

plots = openmc.Plots([plot1, plot2, plot3, plot4, plot5, plot6])
plots.export_to_xml()

# Reaction rates
cell_filter = openmc.CellFilter(fuel_cells[0])
reaction_tally = openmc.Tally(1)
reaction_tally.filters = [cell_filter]
reaction_tally.nuclides = ['U235']
reaction_tally.scores = ['total', 'fission', 'absorption', '(n,gamma)']

# Fuel flux tally
fuel_tally = openmc.Tally()
fuel_tally.filters = [openmc.DistribcellFilter(fuel_cells[0])]
fuel_tally.scores = ['flux']

# Create equal-lethargy energies to put in filter
energies = np.logspace(np.log10(1e-5), np.log10(20.0e6), 501)
e_filter = openmc.EnergyFilter(energies)
# Create tally with energy filter
flux_energy_tally = openmc.Tally()
flux_energy_tally.filters = [e_filter]
flux_energy_tally.scores = ['flux']

openmc.mgxs.GROUP_STRUCTURES.keys()
# Create energy filter using SHEM-361 group structure
energies_shem = openmc.mgxs.GROUP_STRUCTURES['SHEM-361']
shem_filter = openmc.EnergyFilter(openmc.mgxs.GROUP_STRUCTURES['SHEM-361'])
tally_shem = openmc.Tally()
tally_shem.filters = [shem_filter]
tally_shem.scores = ['flux']

tallies = openmc.Tallies([fuel_tally, reaction_tally, tally_shem, flux_energy_tally])
tallies.export_to_xml()
