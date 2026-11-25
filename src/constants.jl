# Physical constants (SI units).

const AVOGADRO = 6.02214076e23              # mol⁻¹
const ELEMENTARY_CHARGE = 1.602176634e-19   # C
const BOLTZMANN = 1.380649e-23              # J K⁻¹

"Faraday constant ``F = N_A e`` (C mol⁻¹)."
const FARADAY = AVOGADRO * ELEMENTARY_CHARGE

"Molar gas constant ``R = N_A k_B`` (J mol⁻¹ K⁻¹)."
const GAS_CONSTANT = AVOGADRO * BOLTZMANN

@doc raw"Vacuum permittivity ``\varepsilon_0`` (F m⁻¹), CODATA 2018 recommended value."
const VACUUM_PERMITTIVITY = 8.8541878128e-12

"""
Relative permittivity of water at 25 °C (78.30), measured by Malmberg and Maryott (1956),
J. Res. Natl. Bur. Stand. 56(1). The temperature dependence is not implemented.
"""
const WATER_RELATIVE_PERMITTIVITY = 78.30
