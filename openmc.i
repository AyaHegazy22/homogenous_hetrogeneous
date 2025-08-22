pi = 3.14159265359
ave_T =  1001.02                           # average temperature (K)

[Mesh]
  type = FileMesh
  file = mesh/multi_layer_equal_volumes_n=8.e
[]


[AuxVariables]
  [material_id]
    family = MONOMIAL
    order = CONSTANT
  []
  [cell_temperature]
    family = MONOMIAL
    order = CONSTANT
  []
  [cell_density]
    family = MONOMIAL
    order = CONSTANT
  []
    [z]
    family = MONOMIAL
    order = CONSTANT
  []
[]



[AuxKernels]
  [material_id]
    type = CellMaterialIDAux
    variable = material_id
  []
  [cell_temperature]
    type = CellTemperatureAux
    variable = cell_temperature
  []
  [cell_density]
    type = CellDensityAux
    variable = cell_density
  []
  [z]
    type = ParsedAux
    variable = z
    use_xyzt = true
    function = 'z'
  []
[]



[ICs]
  [matrix_temp]
    type = ConstantIC
    variable = matrix_temp
    value = ${ave_T}
  []
 [heat_source]
    type = ConstantIC
    variable = heat_source
    value = 157882019.311
    block = '7'
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  verbose = true
  power = 580.5355  
  scaling = 100.0
  temperature_variables = ' matrix_temp ; matrix_temp'
  temperature_blocks = '7 ; 8'
  lowest_cell_level = 1
  check_equal_mapped_tally_volumes = true
  initial_properties = xml
  relaxation = robbins_monro
  particles = 10000
  max_batches = 4000

    [Tallies]
    [cell_tally]
      type = CellTally
      block = '7'
      name = heat_source
      output = 'unrelaxed_tally_std_dev'
      check_equal_mapped_tally_volumes = true
      trigger = rel_err
      trigger_threshold = 1e-3
    []
  []
[]


[Postprocessors]
  [heat_source]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    execute_on = 'transfer initial timestep_end'
    block = '7'
  []
  [max_tally_rel_err]
    type = TallyRelativeError
    value_type = max
  []
  [k]
    type = KEigenvalue
  []
  #[k_std_dev]
  #  type = KStandardDeviation
  #[]
  [min_power]
    type = ElementExtremeValue
    variable = heat_source
    value_type = min
  []
  [max_power]
    type = ElementExtremeValue
    variable = heat_source
    value_type = max
  []
  [z_max_power]
    type = ElementExtremeValue
    proxy_variable = heat_source
    variable = z
  []
[]


[MultiApps]
  [matrix]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'matrix.i'
    execute_on = timestep_begin
  []
[]



[Transfers]
  [source_to_matrix]
    type =  MultiAppGeneralFieldShapeEvaluationTransfer
    source_variable = heat_source
    variable = power
    to_multi_app = matrix
    from_postprocessors_to_be_preserved = heat_source
    to_postprocessors_to_be_preserved = power
  []
  [temp_from_matrix]
    type = MultiAppGeometricInterpolationTransfer
    source_variable = matrix_temp
    from_multi_app = matrix
    variable = matrix_temp
  []
[]

[UserObjects]
  [heat_source_avg]
    type = NearestRadiusLayeredAverage
    num_layers = 1
    points = '0.00035 0 0
              0.00135 0 0
              0.00235 0 0
              0.00335 0 0
              0.00435 0 0
              0.00535 0 0
              0.00635 0 0'
    block = '7'
    direction = z
    variable = heat_source
  []
[] 

[VectorPostprocessors]
  [spatial_from_heat_source]
    type = SpatialUserObjectVectorPostprocessor
    userobject = heat_source_avg
  []
[]


[Executioner]
  type = Transient
  #num_steps = 10
  steady_state_detection = true
  check_aux = true
  steady_state_tolerance = 1e-4
  dt = 10
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
