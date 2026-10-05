"""Single source of numerical inputs for marcus_repro."""
import math

# Charged-lepton pole masses (MeV)
M_E = 0.51099895069
M_MU = 105.6583755
M_TAU = 1776.930

# Quark central values used by the manuscript census (MeV)
R_M = 1.01750
M_U_2GEV = 2.16
M_D_2GEV = 4.70
M_S_2GEV = 92.9
M_C_SELF = 1272.9
M_B_SELF = 4186.0
M_T_DIRECT = 172600.0

M_U = M_U_2GEV * R_M
M_D = M_D_2GEV * R_M
M_S = M_S_2GEV * R_M

NAMES = ["e","mu","tau","u","d","s","c","b","t"]
PDG = [M_E,M_MU,M_TAU,M_U,M_D,M_S,M_C_SELF,M_B_SELF,M_T_DIRECT]

MU_STAR = M_E + M_MU + M_TAU
ALPHA_K = math.sqrt(1.5) - 1.0
G = math.sqrt(1.5) * MU_STAR
P = 2.0 * M_E

# Procedural generators
LEP_LO, LEP_HI = M_E, M_TAU
Q_LO, Q_HI = M_U, M_T_DIRECT
KOIDE_EPS = 2.2e-6

# Observed residual thresholds. These are generated from one declared spectrum,
# except STRANGE_ANCHOR_THRESHOLD which preserves the manuscript's archived
# unrounded 0.6216% value explicitly rather than silently deriving 0.620...%
STRANGE_ANCHOR_THRESHOLD = 0.006216
