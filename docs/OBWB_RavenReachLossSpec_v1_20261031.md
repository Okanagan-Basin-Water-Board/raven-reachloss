# Native losing-reach process for Raven: specification and implementation map (phase A)

**v1, 2026-10-31 (Pacific).** Okanagan Basin Water Board, Water Science and Stewardship. Decisions D179 and D181 in the Mission Creek project decision log. This note is for the Raven developers. It states the requirements, maps each one to the TEST implementation in `patches/OBWB_RavenReachLoss_v1_20261031.diff`, and lists the tests that check it.

## 1. Problem

Raven v4.12 has no loss term in a river reach. Losing reaches are modelled today with `:FlowDiversionLookupTable SBID -1`, which:
- removes the water from the model, although it physically recharges an aquifer;
- has no stage or head form, so it cannot tell a disconnected reach from a connected one;
- in v4.12, prints a stale warning that the watershed balance does not count the loss. The loss is in fact counted in Cum. Outflow (`Model.cpp:2452-2462`; Mission Creek finding F27).

On Mission Creek, the fan groundwater model shows the reaches above WSC 08NM116 are disconnected in every Aug-Sep month, with the aquifer 0.9 to 5.4 m below the bed bottom, while the reaches below are connected (F215). Observed heads agree.

## 2. Requirements and where each is met

| ID | Requirement | Implementation (patch) | Test (`mission_creek/reachloss_mc.py`) |
|---|---|---|---|
| R1 | Mass conserved; lost water has a named destination booked in the watershed balance | HRU and subbasin destinations are netted out of Cum. Outflow (`Model.cpp`); the return store is counted in channel storage (`GetTotalChannelStorage`) | A3: MB error of D2 and D3 within 1 mm of D1 |
| R2 | Per subbasin reach, optional Julian-day window | `:ReachLoss SBID FORM start end`; uses the same `IsInDateRange` as diversions | D1 reproduces D32's window (A2) |
| R3 | Disconnected form, independent of aquifer head | `DISCONNECTED`: Kbed (h + d)/d W length, hydraulic geometry W = aW Q^bW, h = aH Q^bH | A7: law check; D5a equals D4 exactly |
| R3a | Empirical form that reproduces existing practice | `EMPIRICAL_Q`: min(a Q^b, cap) | A2: D1 equals the D32 lookup run |
| R4 | Connected form with head difference, continuous with R3 | `CONNECTED`: Kbed ((z_bed + h) - max(h_gw, z_bed - d))/d W length | A7: D5a = D4 (switch continuity); D5b shows no flow against the gradient |
| R5 | Loss never more than the available water | `min(L, Q)` in `CalcReachLossRate`; actual loss scaled by the supply-limited diversion in `UpdateOutflows` | A4 |
| R6 | Destinations: sink, HRU storage, delayed return | `SINK`; `HRUGROUP name SV` (area-weighted over enabled HRUs, after routing); `SUBBASIN SBID k` (linear store, release to `aRouted`) | A3, A6 |
| R7 | Per-reach output each step | `ReachLoss.csv`: date, SBID, form, destination, Qpre (flow the law used), Qact, vol_step, store, gain clips | A5: volume two ways |
| R8 | Off by default, byte-identical otherwise | No new state, file or arithmetic unless `:ReachLoss` is present (the added storage term is +0.0) | A1: D0 against the reference build, every output file |

## 3. Design choices to discuss with the developers

1. **Built on the diversion machinery.** A reach loss is a `diversion` with `loss_form > 0`, so it uses the existing removal in `UpdateOutflows` and the trapezoid booking (`GetDiversionLosses`). A separate reach-process class may be cleaner upstream.
2. **Head source for `CONNECTED`.** The prototype takes a constant h_gw. Useful extensions are a time series, the level of a Raven storage (for example a `GROUNDWATER` or `AQUIFER` store of the destination HRUs), or the `GWRiverConnection` path used in coupled mode.
3. **Gaining reaches.** The prototype clips gains to zero and counts them. A gain needs a source store with finite volume (requirement P4 of the method).
4. **HRU deposit timing.** The volume is added to `aPhinew` after routing, so the soil processes see it the next step. An alternative is to apply it as a forcing-like flux inside the process solve.
5. **Hydraulic geometry.** W and h use power laws in Q (Leopold and Maddock, 1953). Reaches with a channel profile could use Raven's own stage and width from `CChannelXSect`.
6. **Destination physics matters more than the loss law (Mission Creek result F230).** Sending the loss to the deep soil of the losing reach's HRUs, or to a 56-day return store to the outlet, returned 100% of it to the creek within weeks to months. That cut the kokanee-window low-flow days from 174 to 70 to 73. The fan groundwater model returns only 36% to the creek over 8 years. A useful upstream destination is a slow aquifer store whose discharge can be split between the channel, an outlet such as a lake, and an external sink, each with its own time constant.
7. **The unimplemented `:IrrigationDestination`.** It is parsed but never applied (`DemandOptimization.cpp:2510`, "todo: add irrigation support"). Implementing it would give demand return flows a land destination with the same booking as R6.

## 4. Known limits of the TEST build

- macOS arm64 build only. A Linux container build for the cloud gate is pending.
- `ReachLoss.csv` uses a static stream, so one file per process: ensemble or multi-run modes are not handled.
- No transport (constituent) routing of the lost water. Heat and mass tracers ignore it.
- `:ReachLoss` is parsed from the `.rvt` only.
- The HRU destination does not check the receiving store's maximum capacity. Overfilling is left to Raven's soil processes.

## References

- Leopold, L. B., & Maddock, T. (1953). *The hydraulic geometry of stream channels and some physiographic implications* (Professional Paper 252). U.S. Geological Survey. [training data, High confidence]
- Brunner, P., Cook, P. G., & Simmons, C. T. (2009). Hydrogeologic controls on disconnection between surface water and groundwater. *Water Resources Research, 45*, W01422.
- Raven v4.12 source, commit 81011f0: `SubBasin.cpp`, `Solvers.cpp`, `Model.cpp`, `ParseTimeSeriesFile.cpp`, `DemandOptimization.cpp:2510`.
- OBWB Mission Creek project decision log entries F27, F215, D179, D181 and the method note `docs/OBWB_LosingReachMethod_v1_20261031.md` (private repository; available on request).
