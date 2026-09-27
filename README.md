# LLCC Probabilities

Hourly-updated chance of violating each NASA-STD-4010B Lightning Launch Commit
Criteria rule, and any of them, for LC-39A or SLC-40, in 3-hour windows out 7
days, plus any window you enter (including a near-instantaneous one).

This is a planning aid built from model guidance, not an LLCC evaluation.

## How each rule is estimated

| Rule | Evidence |
|---|---|
| Lightning 4.1.1 | Hi-res lightning field within 10 nmi of the flight path; NBM thunderstorm probability |
| Cumulus 4.1.3 | Reflectivity at the -10 C level within 5 nmi, or echo top above the -20 C height within 10 nmi |
| Anvil 4.1.4-4.1.5 (rough) | Deep convection within 20 nmi with high cloud within 10 nmi |
| Debris 4.1.6 (rough) | Mid/high cloud, no deep convection now, deep convection within 10 nmi in the previous 3 h |
| Disturbed weather 4.1.7 (rough) | 30 dBZ or more within 5 nmi under an overcast mid deck |
| Thick cloud layers 4.1.8 (rough) | Mid cloud near overcast with echo within 5 nmi |
| Surface electric fields, smoke, triboelectrification | Not forecast (see page notes) |

Direct evidence comes from HRRR, the HREF member models, NAM 3 km and NBM.
Where none covers a window (mostly beyond ~2.5 days) the global ensembles
supply stand-ins (deep convection = precip with CAPE >= 1000 J/kg), shown hatched.

The flight path is a line from the pad along `azimuth_deg` for `length_nmi`
(`config.yaml`). The -20 C height (`h20_m`) is a seasonal setting; update it
as the season changes.

## Set up

Same as the other boards: new repo, push with git, Pages from `/docs`,
Actions write permission, run the workflow once and read the probe log. The
probe lists which of the needed fields (LTNG, REFD at 263 K, RETOP, HCDC,
MCDC, REFC) each hi-res model provides; rules that need a missing field skip
that model.
