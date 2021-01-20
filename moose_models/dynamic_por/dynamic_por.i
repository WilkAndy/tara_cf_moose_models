# With dynamic porosity and dependence on volumetric strain rate
[Mesh]
  type = GeneratedMesh
  dim = 2
# nx=50, bias_x=1.2
# nx=100, bias_x=1.1
# nx=200, bias_x=1.05
# nx=400, bias_x=1.02
# nx=1000, bias_x=1.01
# nx=2000, bias_x=1.003
  nx = 2000
  bias_x = 1.003
  xmin = 0.1
  xmax = 5000
  ny = 1
  ymin = 0
  ymax = 11
[]

[Problem]
  coord_type = RZ
[]
  
[GlobalParams]
  displacements = 'disp_r disp_z'
  PorousFlowDictator = dictator
  gravity = '0 0 0'
  biot_coefficient = 1.0
[]

[Variables]
  [./pwater]
    initial_condition = 18.3e6
  [../]
  [./sgas]
    initial_condition = 0.0
  [../]
  [./temp]
    initial_condition = 358
  [../]
  [./disp_r]
  [../]
[]

[AuxVariables]
  [./rate]
  [../]
  [./disp_z]
  [../]
  [./massfrac_ph0_sp0]
    initial_condition = 1 # all H20 in phase=0
  [../]
  [./massfrac_ph1_sp0]
    initial_condition = 0 # no H2O in phase=1
  [../]
  [./pgas]
    family = MONOMIAL
    order = FIRST
  [../]
  [./swater]
    family = MONOMIAL
    order = FIRST
  [../]
  [./stress_rr]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_tt]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  [../]
  [./flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  [../]
  [./vol_strain_rate_water]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  [../]
  [./mass_co2_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    use_displaced_mesh = false
    variable = sgas
  [../]
  [./flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    use_displaced_mesh = false
    variable = sgas
  [../]
  [./vol_strain_rate_co2]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 1
    use_displaced_mesh = false
    variable = sgas
  [../]
  [./energy_dot]
    type = PorousFlowEnergyTimeDerivative
    use_displaced_mesh = false
    variable = temp
  [../]
  [./advection]
    type = PorousFlowHeatAdvection
    use_displaced_mesh = false
    variable = temp
  [../]
  [./conduction]
    type = PorousFlowExponentialDecay
    use_displaced_mesh = false
    variable = temp
    reference = 358
    rate = rate
  [../]
  [./vol_strain_rate_heat]
    type = PorousFlowHeatVolumetricExpansion
    use_displaced_mesh = false
    variable = temp
  [../]
  [./grad_stress_r]
    type = StressDivergenceRZTensors
    temperature = temp
    thermal_eigenstrain_name = thermal_contribution
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  [../]
  [./poro_r]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  [../]
[]

[AuxKernels]
  [./rate]
    type = FunctionAux
    variable = rate
    execute_on = timestep_begin
    function = decay_rate
  [../]
  [./pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas
  [../]
  [./swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  [../]
  [./stress_rr]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_rr
    index_i = 0
    index_j = 0
  [../]
  [./stress_tt]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_tt
    index_i = 2
    index_j = 2
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 1
    index_j = 1
  [../]
[]

[Functions]
  [./decay_rate]
# Eqn(26) of the first paper of Tara
# Ka * (rho C)_a = 10056886.914
# h = 11
    type = ParsedFunction
    value = 'sqrt(10056886.914/t)/11.0'
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'temp pwater sgas disp_r'
    number_fluid_phases = 2
    number_fluid_components = 2
  [../]
  [./pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  [../]
[]

[Modules]
  [./FluidProperties]
    [./water]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 970.0
      viscosity = 0.3394e-3
      cv = 4149.0
      cp = 4149.0
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    [../]
    [./co2]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 516.48
      viscosity = 0.0393e-3
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
    temperature = temp
    at_nodes = true
  [../]
  [./temperature_qp]
    type = PorousFlowTemperature
    temperature = temp
  [../]
  [./ppss]
    type = PorousFlow2PhasePS
    at_nodes = true
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc
  [../]
  [./ppss_qp]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  [../]
  [./water_at_nodes]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
    at_nodes = true
  [../]
  [./water_at_qp]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  [../]
  [./gas_at_nodes]
    type = PorousFlowSingleComponentFluid
    fp = co2
    phase = 1
    at_nodes = true
  [../]
  [./gas_at_qp]
    type = PorousFlowSingleComponentFluid
    fp = co2
    phase = 1
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    at_nodes = true
    include_old = true
    material_property = PorousFlow_fluid_phase_density_nodal
  [../]
  [./dens_all_at_quadpoints]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_nodes = false
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_viscosity_nodal
  [../]
  [./energy_all]
    type = PorousFlowJoiner
    include_old = true
    at_nodes = true
    material_property = PorousFlow_fluid_phase_internal_energy_nodal
  [../]
  [./enthalpy_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_fluid_phase_enthalpy_nodal
  [../]
  [./porosity_reservoir]
    type = PorousFlowPorosity
    porosity_zero = 0.204284486 # ensures porosity=0.2 at T=358
    thermal_expansion_coeff = 15E-6 # 3*5E-6 because this is the volumetric one
    solid_bulk = 8E9 # unimportant since biot=1
    at_nodes = true
  [../]
  [./permeability_reservoir]
    type = PorousFlowPermeabilityConst
    permeability = '2e-12 0 0  0 0 0  0 0 0'
  [../]
  [./relperm_liquid]
    type = PorousFlowRelativePermeabilityCorey
    at_nodes = true
    n = 4
    phase = 0
    s_res = 0.200
    sum_s_res = 0.405
  [../]
  [./relperm_gas]
    type = PorousFlowRelativePermeabilityBC
    at_nodes = true
    phase = 1
    s_res = 0.205
    sum_s_res = 0.405
    nw_phase = true
    lambda = 2
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_relative_permeability_nodal
  [../]
  [./thermal_conductivity_reservoir]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '0 0 0  0 1.320 0  0 0 0'
    wet_thermal_conductivity = '0 0 0  0 3.083 0  0 0 0'
  [../]
  [./internal_energy_reservoir]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1100
    density = 2350.0
  [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 6.0E9
    poissons_ratio = 0.2
  [../]
  [./strain]
    type = ComputeAxisymmetricRZSmallStrain
    eigenstrain_names = 'thermal_contribution ini_stress'
  [../]
  [./ini_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '-12.8E6 0 0  0 -51.3E6 0  0 0 -12.8E6'
    eigenstrain_name = ini_stress
  [../]
  [./thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temp
    stress_free_temperature = 358
    thermal_expansion_coeff = 5E-6
    eigenstrain_name = thermal_contribution
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    at_nodes = true
  [../]
  [./eff_fluid_pressure_qp]
    type = PorousFlowEffectiveFluidPressure
  [../]
  [./vol_strain]
    type = PorousFlowVolumetricStrain
  [../]
[]

[BCs]
  [./outer_pressure_fixed]
    type = DirichletBC
    boundary = right
    value = 18.3e6
    variable = pwater
  [../]
  [./outer_saturation_fixed]
    type = DirichletBC
    boundary = right
    value = 0.0
    variable = sgas
  [../]
  [./outer_temp_fixed]
    type = DirichletBC
    boundary = right
    value = 358
    variable = temp
  [../]
  [./co2_injection]
    type = PorousFlowSink
    boundary = left
    variable = sgas
    use_mobility = false
    use_relperm = false
    fluid_phase = 1
    flux_function = 'min(t/100.0,1)*(-2.294001475)' # 5.0E5 T/year = 15.855 kg/s, over area of 2Pi*0.1*11
  [../]
  [./cold_co2]
    type = PresetBC
    boundary = left
    variable = temp
    value = 294
  [../]
  [./cavity_pressure_x]
    type = Pressure
    boundary = left # probably better than injection_area
    variable = disp_r
    component = 0
    postprocessor = p_bh # note, this lags
    use_displaced_mesh = false
  [../]
  [./fixed_outer_r]
    type = PresetBC
    variable = disp_r
    value = 0
    boundary = right
  [../]
[]

[Postprocessors]
  [./p_bh]
    type = PointValue
    variable = pwater
    point = '0.1 0 0'
    execute_on = timestep_begin
    use_displaced_mesh = false
  [../]
[]

[VectorPostprocessors]
  [./ptsuss]
    type = LineValueSampler
    use_displaced_mesh = false
    start_point = '0.1 0 0'
    end_point = '5000 0 0'
    sort_by = x
    num_points = 50000
    outputs = csv
    variable = 'pwater temp sgas disp_r stress_rr stress_tt'
  [../]
[]
  
[Preconditioning]
  active = 'mumps'
  [./smp]
    type = SMP
    full = true
    #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E2       1E-5        500'
  [../]
  [./mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-5       1E2       50'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 1.5768e8
  #dtmax = 1e6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.1
  [../]
[]

[Outputs]
  print_linear_residuals = false
  sync_times = '3600 86400 2.592E6 1.5768E8'
  print_perf_log = true
  [./exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
    sync_only = true
  [../]
[]
