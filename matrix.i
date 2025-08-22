inlet_T = 598.0                          # inlet fluid temperature (K)
buffer_k = 0.5                           # buffer thermal conductivity (W/m/K)
PyC_k = 4.0                              # PyC thermal conductivity (W/m/K)
SiC_k = 13.9                             # SiC thermal conductivity (W/m/K)
kernel_k = 3.5                           # fissil kernel thermal conductivity (W/m/K)
matrix_k= 15.0                          # graphite matrix thermal conductivity (W/m/K)
kernel_radius =  290.0e-6                # fissile kernel outer radius (m)
buffer_radius = 370.0e-6                # buffer outer radius (m)
iPyC_radius = 410.0e-6                  # inner PyC outer radius (m)
SiC_radius = 445.0e-6                   # SiC outer radius (m)
oPyC_radius = 485.0e-6                  # outer PyC outer radius (m)
height = 0.025                            # height (m)
triso_pf = 0.4

ave_T = 1001.02                          # average temperature (K)


# compute the volume fraction of each TRISO layer in a TRISO particle
# for use in computing average thermophysical properties
kernel_fraction = ${fparse kernel_radius^3 / oPyC_radius^3}
buffer_fraction = ${fparse (buffer_radius^3 - kernel_radius^3) / oPyC_radius^3}
ipyc_fraction = ${fparse (iPyC_radius^3 - buffer_radius^3) / oPyC_radius^3}
sic_fraction = ${fparse (SiC_radius^3 - iPyC_radius^3) / oPyC_radius^3}
opyc_fraction = ${fparse (oPyC_radius^3 - SiC_radius^3) / oPyC_radius^3}


[Mesh]
  type = FileMesh
  file = mesh/multi_layer_equal_volumes_n=8.e
[]


[Variables]
  [matrix_temp]
    initial_condition = ${ave_T}
  []
[]

[Kernels]
  [diffusion]
    type = HeatConduction
    variable = matrix_temp
    block = '7 8'
  []
  [source]
    type = CoupledForce
    variable = matrix_temp
    v = power
    block = '7'   #fuel kernel
  []
[]


[BCs] # Make it neumann on top and bottom and sides, and Dirichlet on the parts touching the coolant
  [top_bottom]
    type = NeumannBC
    variable =  matrix_temp
    boundary = '120 30 40 60 110 70 90'
  []
  [surface]
    type = DirichletBC
    variable =  matrix_temp
    boundary = '100 50 80'
    value =  ${ave_T}
  []
[]

[Functions]
 [k_graphite]
    type = ParsedFunction
    expression = '${matrix_k}'
  []
  [k_TRISO]
    type = ParsedFunction
    expression = '${kernel_fraction} * ${kernel_k} + ${buffer_fraction} * ${buffer_k} + ${fparse ipyc_fraction + opyc_fraction} * ${PyC_k} + ${sic_fraction} * ${SiC_k}'
  []
  [k_compacts]
    type = ParsedFunction
    expression = '${triso_pf} * k_TRISO + ${fparse 1.0 - triso_pf} * k_graphite'
    symbol_names = 'k_TRISO k_graphite'
    symbol_values = 'k_TRISO k_graphite'
  []
[]

[Materials]
  [graphite]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = k_graphite
    temp = matrix_temp
    block = '8'
  []
  [compacts]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = k_compacts
    temp = matrix_temp
    block = '7'
  []
[]

[AuxVariables]
  [power]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 157882019.311
  []
   [temp_average]
    order = CONSTANT
    family = MONOMIAL
  [] 
[]

[AuxKernels]
  [temp_average]
    type = SpatialUserObjectAux
    variable = temp_average
    execute_on = timestep_end
    user_object = Temphistogram
  []
[]

[UserObjects]
  [Temphistogram]
    type = NearestPointAverage
    points_file = points.txt
    variable = matrix_temp
  []
[]


[Postprocessors]
  [power] # evaluate the total power for normalization
    type = ElementIntegralVariablePostprocessor
    variable = power
    execute_on = 'transfer'
    block = '7'
  []
  [max_matrix_T]
    type = ElementExtremeValue
    variable = matrix_temp
    value_type = max
    block = '7 8'
  []
  [avg_matrix_T]
    type = ElementAverageValue
    variable = matrix_temp
    block = '7 8'
  []
[]


[Executioner]
  type = Transient
  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-5
  petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  sub_cycling = true
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
