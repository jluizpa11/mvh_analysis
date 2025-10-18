[paths]
filepath = "examples/example_titration.csv"

[experiment]
wavelength_nm = 506
conc_DNA_uM = 1150.0
conc_dye_uM = 44.0
initial_volume_uL = 500.0

[fit]
p0 = [1e6, 2.0]
bounds_lower = [1e3, 0.5]
bounds_upper = [1e8, 10.0]
function = "mvh_rhs"

[template]
template = "examples/template.csv"
