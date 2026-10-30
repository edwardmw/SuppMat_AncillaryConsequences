# Replication Files

The files required for replicating the main results and the sensitivity analyses in the Supplementary Information (SI) of "Ancillary Consequences of Targeted Policy Interventions in the Presence of Disease-Based Poverty Traps."

The EEL-code.gms file was produced using GAMS v24.1.3 and should work with later versions.

## Files needed

- `EEL-code.gms` (the model)
- `EEL-data.gdx` (the input data, must be in the same folder as EEL-code.gms)
- An `output` subfolder in the same directory, created before running (GAMS will not create it automatically)
- `output_SI_Lag` and `output_SI_Migration` subfolders, created before running, only needed for the sensitivity analyses described below

## Running a scenario

Scenarios are chosen with command-line flags passed to GAMS when launching the run, rather than by editing the code. The base command is

```bash
gams EEL-code.gms --RoWtrade=1 --MDAshock=1 --FMPshock=0 --TFPshock=0
```

If `gams` is not recognized as a command, you may need to add your GAMS installation's directory to your system PATH, or call the full path to the `gams` executable directly.

Each flag is either 0 (off) or 1 (on).

- `RoWtrade` chooses the trade regime for fish. Set to 1 for the Rest-of-World (RoW) trade scenario, the paper's preferred specification. Set to 0 for the Local trade scenario. If this flag is omitted, GAMS defaults to 0 (Local trade), so RoWtrade=1 must be passed explicitly to reproduce the RoW results reported in the main text.
- `MDAshock`, `FMPshock`, and `TFPshock` turn the mass drug administration, fishery management, and total factor productivity interventions on or off. Any combination can be active at once, which is how the paper's concurrent-policy results are produced.

To run a different scenario, change the value of any flag to 0 or 1 and rerun the command. For example, the seven single-domain and concurrent RoW policy combinations reported in Table 1 are produced by

```bash
gams EEL-code.gms --RoWtrade=1 --MDAshock=1 --FMPshock=0 --TFPshock=0
gams EEL-code.gms --RoWtrade=1 --MDAshock=0 --FMPshock=1 --TFPshock=0
gams EEL-code.gms --RoWtrade=1 --MDAshock=0 --FMPshock=0 --TFPshock=1
gams EEL-code.gms --RoWtrade=1 --MDAshock=1 --FMPshock=1 --TFPshock=0
gams EEL-code.gms --RoWtrade=1 --MDAshock=1 --FMPshock=0 --TFPshock=1
gams EEL-code.gms --RoWtrade=1 --MDAshock=0 --FMPshock=1 --TFPshock=1
gams EEL-code.gms --RoWtrade=1 --MDAshock=1 --FMPshock=1 --TFPshock=1
```

Table 2's Local-trade results use the same seven combinations with `--RoWtrade=0`.

## Output

Each run writes a single GDX file to the `output` folder, named after the active flags. For example, `--RoWtrade=1 --MDAshock=1 --FMPshock=0 --TFPshock=0` produces `output/RoW_MDA.gdx`, and `--RoWtrade=0 --FMPshock=1 --TFPshock=1` produces `output/Loc_FMP_TFP.gdx`. Running with every shock flag at 0 produces the baseline, with no policy applied.

## Sensitivity analyses in the Supplementary Information

The sensitivity analyses are switched on with additional flags. All of them are off by default, so the commands above still reproduce the main results exactly. They are combined with the shock and trade flags described above.

### Labour supply elasticity (SI Appendix A6)

`SI_LabElas1=1` sets the labour supply elasticity to 1 instead of 100. `SI_LabElasINF=1` makes the elasticity decline with aggregate infection. For example,

```bash
gams EEL-code.gms --RoWtrade=1 --FMPshock=1 --SI_LabElas1=1
```

produces `output/RoW_FMP_LabElas1.gdx`. The second specification adds the suffix `_LabElasINF` instead.

### Lagged infection burden (SI Appendix A7)

`SI_Lag=1` makes the infection burden that enters effective labour adjust only partially toward the equilibrium prevalence each year. `lagLambda` is the share of the gap closed each year and `lagtag` is the label added to the output file name. A value of 1 for `lagLambda` means no lag and reproduces the main results. The lag is not applied when `MDAshock=1`, because MDA clears adult worms directly. For example,

```bash
gams EEL-code.gms --RoWtrade=1 --FMPshock=1 --SI_Lag=1 --lagLambda=0.25 --lagtag=025
```

produces `output_SI_Lag/RoW_FMP_Lag025.gdx`. The paper uses values of 0.5, 0.25 and 0.1 (tags 050, 025 and 010) for the FMP and TFP interventions under both trade regimes.

### Net migration (SI Appendix A8)

`SI_Mig=1` lets the population of each household type respond to its own income per capita. `migEta` is the annual elasticity of population with respect to income per capita and `migtag` is the label added to the output file name. `migAll=1` applies the rule to all four household types, as in the paper, while the default of 0 applies it to the two fishing household types only. Migrants bring labour only. For example,

```bash
gams EEL-code.gms --RoWtrade=1 --FMPshock=1 --SI_Mig=1 --migAll=1 --migEta=0.05 --migtag=a050
```

produces `output_SI_Migration/RoW_FMP_Miga050.gdx`. The paper uses values of 0.025 and 0.05 (tags a025 and a050) for the MDA, FMP and TFP interventions under both trade regimes. Combining these flags with `--SI_LabElas1=1` runs the migration analysis at a labour supply elasticity of 1, and the file names then carry the extra label `_LabElas1`.

Three further flags control how population enters the epidemiological side of the model. They are used to decompose the effect of migration on the basic reproduction number.

- `migExpPC` (default 0). When set to 1, each household's fishing labour time is divided by its population index inside the exposure ratio, which gives exposure per household member. The paper's results use 0, so exposure is the household's fishing labour time relative to its baseline value.
- `migShares` (default 1). Updates the population shares that weight household types each year. Setting 0 holds them at their baseline values.
- `migBG` (default 1). Lets population enter the output per capita that determines the transmission and recovery rates, beta and gamma. Setting 0 holds population at its baseline value there.

### Checks

Setting `lagLambda=1` or `migEta=0` reproduces the main results. This can be confirmed by comparing the resulting GDX file with the corresponding file in `output` using the GAMS utility `gdxdiff`.
