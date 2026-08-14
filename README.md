# Replication Files

The files required for replicating the main results presented in "Ancillary Consequences of Targeted Policy Interventions in the Presence of Disease-Based Poverty Traps."

The EEL-code.gms file was produced using GAMS v24.1.3 and should work with later versions.

## Files needed

- `EEL-code.gms` (the model)
- `EEL-data.gdx` (the input data, must be in the same folder as EEL-code.gms)
- An `output` subfolder in the same directory, created before running (GAMS will not create it automatically)

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
