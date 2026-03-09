# Energy Storage State-of-Charge Management in Rolling Horizon Markets
### A Pure Pyomo Implementation on the RTS-GMLC Test System

---

## Overview

This project implements a rolling horizon electricity market simulation with explicit State-of-Charge (SoC) management for grid-scale energy storage resources (ESRs). The core idea is taken from market operations research: in a real electricity market, grid operators cannot see the full day ahead — they solve short-horizon optimization problems repeatedly, updating forecasts as the day unfolds. This creates a fundamental tension for energy-limited resources like batteries: discharge now to handle current net load variability, or preserve energy for expected future needs?

This implementation explores three SoC management strategies — **Greedy**, **Conservative**, and **Balanced** — and compares them against a **Full Horizon** perfect foresight benchmark. The key finding being validated is that a moderate penalty on SoC target violations produces reliability outcomes equivalent to perfect foresight, while naive strategies (too aggressive or too conservative) worsen system reliability.

The implementation is entirely in **Python + Pyomo** using **free open-source solvers (HiGHS/CBC)** and the publicly available **RTS-GMLC test system**, making it fully reproducible without commercial software or proprietary data.

---

## Motivation

Most academic studies on storage in electricity markets assume the optimizer has perfect knowledge of the entire day's net load. Real markets operate very differently — operators solve a rolling horizon with constantly updating short-term forecasts. This gap between academic models and operational reality has important implications for how storage should be dispatched.

The central challenge is SoC management under uncertainty:

```
Net Load = Total Demand - Wind - Solar - Nuclear - Hydro

Operators don't know full day net load because:
→ Demand forecast errors accumulate over the day
→ Wind output is uncertain beyond a few hours
→ Solar output depends on cloud cover — unpredictable at short timescales
→ All three errors compound into net load forecast error
```

Two naive strategies both fail:

- **Greedy** — discharge storage aggressively whenever a net load spike occurs. Storage depletes early in the day and has nothing left for the evening peak. Reliability is actually worse than having no storage at all.
- **Conservative** — hold storage strictly to the day-ahead schedule. Storage sits idle during unexpected mid-day volatility and misses short-term reliability events.

The **Balanced** strategy navigates between these by using a moderate penalty on SoC target violations — allowing deviations when the reliability benefit is large enough, but discouraging unnecessary depletion.

---

## Market Simulation Structure

The simulation replicates the nested structure of North American real-time electricity markets:

```
Day-Ahead (DA)
  Horizon:    24 hours
  Resolution: 1 hour
  Model type: MILP
  Purpose:    Commit slow-start thermal units
              Set hourly SoC schedule for storage
              Produces: u_g,t (commitment), e_DA_s,t (SoC targets)
        ↓
Real-Time Unit Commitment (RTUC)
  Horizon:    2 hours (8 x 15-min intervals)
  Resolution: 15 minutes
  Frequency:  Every 15 minutes throughout the day
  Model type: MILP
  Purpose:    Commit fast-start units
              Re-optimize SoC within 2hr window
              Terminal SoC must stay near e_DA target (soft constraint)
        ↓
Real-Time Economic Dispatch (RTED)
  Horizon:    Single interval
  Resolution: 5 minutes
  Frequency:  Every 5 minutes (288 times per day)
  Model type: LP
  Purpose:    Dispatch all committed units
              Storage charge/discharge decision
              Terminal SoC must stay near e_RTUC target (soft constraint)
              PENALTY ON VIOLATION = strategy selector
        ↓
Real-Time Pricing Run (RTPR)
  Horizon:    Single interval
  Resolution: 5 minutes
  Frequency:  After each RTED
  Model type: LP (binaries relaxed)
  Purpose:    Compute real-time LMPs
              Fast-start startup costs reflected in prices
```

---

## The Three Strategies

All three strategies use **identical model structure and constraints**. The only difference is one parameter — the penalty cost for violating the SoC target in RTED:

```python
strategies = {
    'greedy':        {'penalty': 50},    # $/MWh — cheap to violate
    'balanced':      {'penalty': 220},   # $/MWh — moderate
    'conservative':  {'penalty': 11000}, # $/MWh — extremely expensive
    'full_horizon':  {'perfect_foresight': True}
}
```

### Why One Parameter Controls Everything

The RTED objective is:

```
Minimize:
    generation cost
  + storage variable O&M
  + unserved load penalty
  + penalty_price * SoC_violation_slack    ← this term controls strategy
```

The SoC target is a soft constraint:

```
SoC[s, end_of_interval] + slack[s] >= SoC_target[s]
slack[s] >= 0
cost of slack = slack * penalty_price
```

At each interval the optimizer compares:
- Cost of serving load by discharging storage now
- Cost of violating SoC target

```
If penalty is low ($50):
→ cheap to violate target → discharge freely → Greedy

If penalty is high ($11,000):
→ expensive to violate → hold storage rigidly → Conservative

If penalty is moderate ($220):
→ discharges when reliability benefit exceeds penalty
→ preserves enough SoC for later needs → Balanced
```

---

## Mathematical Formulation

All equations implemented directly from the paper in Pyomo.

### Objective Function

```
min  Σ_g,t c_g,t  +  Σ_s,t C_vom*(p_dis + p_ch)  +  C_uls*ul_t  +  C_sr*s_sr_j,t
```

| Term | Description |
|---|---|
| c_g,t | Thermal generation cost including startup/shutdown |
| C_vom*(p_dis+p_ch) | Storage variable O&M — cost of cycling |
| C_uls * ul_t | Unserved load penalty |
| C_sr * s_sr | Reserve shortage penalty |

### Core Constraints

**Energy Balance (every time period):**
```
Σ p_g,t  +  Σ p_dis_s,t  -  Σ p_ch_s,t  +  ul_t  >=  D_t
```

**Reserve Balance:**
```
Σ r_sr_g,t  +  Σ p_sr_s,t  +  s_sr_t  >=  R_sr
```

**SoC Evolution (battery physics):**
```
e_s,t = e_s,t-1  +  η_ch * p_ch_s,t * τ  -  p_dis_s,t * τ / η_dis
```

**SoC Limits:**
```
E_min_s  <=  e_s,t  <=  E_max_s
```
(Prevents deep cycling — protects battery lifetime)

**Charge/Discharge Power Limits:**
```
P_min * u_ch_s,t  <=  p_ch_s,t  <=  P_max * u_ch_s,t
P_min * u_dis_s,t <=  p_dis_s,t <=  P_max * u_dis_s,t
```

**Mutual Exclusivity (cannot charge and discharge simultaneously):**
```
u_ch_s,t  +  u_dis_s,t  <=  1
```

**Reserve Qualification (must have enough SoC to deliver reserves):**
```
p_sr_s,t * T  <=  (e_s,t - E_min_s) * η_dis
```

**Joint Capacity (energy + reserves cannot exceed power rating):**
```
p_sr_s,t  +  p_dis_s,t  -  p_ch_s,t  <=  P_max_s
```

**DA Terminal SoC:**
```
e_s,t = TCE_s    for t = t_start and t = t_end
```

### Rolling Horizon SoC Constraints

**RTUC — Initial condition:**
```
e_s,t = e_RTED_s,t    (start from actual last RTED SoC)
```

**RTUC — Soft terminal target:**
```
Σ_s e_s,t_end  +  slack_ue  >=  Σ_s e_DA_s,t_end
cost of slack_ue = C_ue * slack_ue   (added to objective)
```

**RTED — Initial condition:**
```
e_s,t = e_RTED_s,t    (start from actual last RTED SoC)
```

**RTED — Soft terminal target (THE STRATEGY KNOB):**
```
Σ_s e_s,t_end  +  slack_ue  >=  Σ_s e_RTUC_s,t_end
cost of slack_ue = penalty_price * slack_ue

penalty_price = 50     → Greedy
penalty_price = 220    → Balanced
penalty_price = 11000  → Conservative
```

---

## Data

### Source: RTS-GMLC Test System

```
github.com/GridMod/RTS-GMLC
```

The RTS-GMLC (Reliability Test System — Grid Modernization Lab Consortium) is a publicly available, modern power systems test case designed specifically for production cost modeling with renewables and storage.

```
Network:      73 buses, 3 areas
Branches:     120 transmission lines
Generators:   158 units (thermal, wind, solar, hydro)
Time series:  8760 hourly + 5-minute resolution
              Load profiles per bus
              Wind generation profiles
              Solar generation profiles
              Reserve requirements
```

### What We Use From RTS-GMLC

```
gen.csv           → generator parameters (pmin, pmax, costs, ramp rates, min up/down)
bus.csv           → bus data
branch.csv        → transmission lines and reactances
timeseries/LOAD_* → hourly load per bus
timeseries/Wind_* → wind generation profiles
timeseries/Solar_*→ solar generation profiles
```

### Storage Units (Defined Manually)

Storage is not natively in RTS-GMLC in the configuration we need, so we define a small fleet scaled to the test system size:

```python
storage_fleet = {
    'BESS_1': {'p_max': 100, 'energy_cap': 400, 'duration': 4},  # 4-hour BESS
    'BESS_2': {'p_max': 100, 'energy_cap': 400, 'duration': 4},
    'BESS_3': {'p_max': 150, 'energy_cap': 900, 'duration': 6},  # 6-hour BESS
    'BESS_4': {'p_max': 150, 'energy_cap': 900, 'duration': 6},
}

# All units:
efficiency     = 0.96       # 96% round-trip
soc_min        = 0.20       # 20% minimum
soc_max        = 0.80       # 80% maximum
vom_cost       = 3.0        # $/MWh variable O&M
```

### Forecast Error Model

Real-time forecasts diverge from the day-ahead forecast throughout the day. We simulate this with an ARIMA-based error model:

```python
# Fit ARIMA to synthetic forecast errors
# Generates three diverging forecast streams:
#   DA forecast    → generated at midnight (highest error)
#   RTUC forecast  → updated every 15min (moderate error)
#   RTED actuals   → ground truth (what actually happens)
```

Simple Gaussian approximation is used initially, with ARIMA as an upgrade path.

---

## Validation Metrics

We validate the implementation by reproducing the key result from the paper's logic: **the Balanced strategy should eliminate unserved energy and match Full Horizon (perfect foresight) outcomes.**

### Primary Metrics

| Metric | Description | Expected Result |
|---|---|---|
| Unserved Energy (GWh) | Load not served due to insufficient supply | Balanced = 0, Greedy > No Storage |
| Reserve Shortage (GWh) | Reserve requirements not met | Balanced = 0 |
| Total System Cost ($) | Fuel + startup + storage VOM | Balanced < Conservative |
| SoC Trajectory | Actual SoC vs DA target vs RTUC target | Balanced tracks DA without depleting |

### Strategy Comparison Table (Target)

| Strategy | Unserved Energy | Reserve Shortage | Behavior |
|---|---|---|---|
| No Storage | > 0 | > 0 | Baseline reference |
| Greedy | > No Storage | > No Storage | Depletes early, fails at peak |
| Conservative | > 0 | ~0 | Holds too long, misses mid-day events |
| Balanced | 0 | 0 | Handles both short-term and peak |
| Full Horizon | 0 | 0 | Perfect foresight benchmark |

### Secondary Metrics

- **LMP time series** — real-time prices from RTPR pricing run
- **Storage dispatch stack** — charge/discharge profile throughout the day
- **SoC divergence** — gap between DA target and actual RTED SoC
- **Penalty activation frequency** — how often SoC slack is non-zero per strategy

---

## Project Structure

```
soc_market_sim/
├── data/
│   ├── rts_gmlc/          ← raw RTS-GMLC CSV files (downloaded from GitHub)
│   └── json/              ← processed EGRET-compatible JSON files
│       ├── DA_data.json
│       └── RT_data.json
│
├── src/
│   ├── input_manager/
│   │   ├── csv_reader.py       ← read RTS-GMLC CSVs
│   │   ├── json_builder.py     ← convert to EGRET schema
│   │   └── forecast_generator.py ← generate diverging forecast streams
│   │
│   ├── models/
│   │   ├── da_model.py         ← DA MILP (24hr, hourly)
│   │   ├── rtuc_model.py       ← RTUC MILP (2hr, 15-min)
│   │   ├── rted_model.py       ← RTED LP (single 5-min interval)
│   │   └── rtpr_model.py       ← RTPR LP (pricing run)
│   │
│   ├── constraints/
│   │   ├── thermal_uc.py       ← standard UC constraints
│   │   ├── storage_soc.py      ← SoC evolution, limits, joint capacity
│   │   └── soc_penalty.py      ← soft SoC target + penalty framework
│   │
│   ├── simulator/
│   │   ├── rolling_horizon.py  ← main simulation loop
│   │   ├── state_manager.py    ← pass SoC state between solves
│   │   └── strategy.py         ← greedy / balanced / conservative config
│   │
│   └── output_manager/
│       ├── result_extractor.py ← pull variables from solved Pyomo model
│       ├── metrics.py          ← compute unserved energy, reserve shortage
│       └── plots.py            ← SoC trajectory, dispatch stack, LMP plots
│
├── notebooks/
│   ├── 01_data_exploration.ipynb    ← explore RTS-GMLC data
│   ├── 02_da_model_validation.ipynb ← validate DA UC on 10-unit system
│   ├── 03_storage_validation.ipynb  ← validate SoC constraints
│   └── 04_strategy_comparison.ipynb ← main results
│
├── tests/
│   ├── test_da_model.py     ← 10-unit system, compare cost to known result
│   ├── test_soc_evolution.py ← verify battery physics
│   └── test_rolling_horizon.py
│
├── requirements.txt
└── README.md
```

---

## Implementation Phases

### Phase 1 — DA Model + UC Validation 
- Build DA MILP in pure Pyomo
- Thermal generators only, no storage
- Validate on classic 10-unit system
- Cost should match Kazarlis 1996 benchmark (~$563,977)
- Then move to RTS-GMLC

### Phase 2 — Add Storage to DA 
- Implement SoC evolution constraint (Eq 4)
- Implement SoC limits (Eq 5)
- Implement power limits and mutual exclusivity (Eq 6-8)
- Implement reserve qualification and joint capacity (Eq 9-10)
- Implement terminal SoC (Eq 11)
- Verify charge/discharge behavior on simple 3-bus example

### Phase 3 — RTED Single Interval
- Build LP for single 5-minute interval
- Commitment fixed as parameter from DA
- Implement soft SoC target (Eq 14-15)
- Test penalty knob manually — does behavior change as expected?

### Phase 4 — RTUC 
- Build MILP for 2-hour, 15-minute horizon
- Fast-start unit commitment
- Implement soft terminal SoC toward DA target (Eq 12-13)
- Connect RTUC output to RTED input

### Phase 5 — Rolling Horizon Loop 
- Build Python control flow: DA → RTUC every 15min → RTED every 5min
- Pass SoC state between solves correctly
- Pass SoC targets: DA → RTUC → RTED
- Add forecast divergence at each step

### Phase 6 — Strategy Comparison + Results 
- Run all four strategies on same representative day
- Compute unserved energy, reserve shortage, total cost
- Plot SoC trajectories (reproduce paper Figure 3 logic)
- Plot net load divergence (reproduce paper Figure 2 logic)
- Build results table (reproduce paper Table 1 logic)

---

## Dependencies

```
python >= 3.9
pyomo >= 6.0
highspy          ← free MIP/LP solver
pandas
numpy
matplotlib
statsmodels      ← for ARIMA forecast error model
egret            ← network data and PTDF (optional, for network layer)
```

Install:
```bash
pip install pyomo highspy pandas numpy matplotlib statsmodels
pip install egret  # optional
```

---

## Key Design Decisions



**Why RTS-GMLC and not a larger system?**
RTS-GMLC provides everything needed — network, generators, full year time series for load, wind, and solar — and runs on a standard laptop. The structural result (Balanced strategy outperforms Greedy and Conservative) does not depend on network size.

**Why HiGHS and not CPLEX?**
Reproducibility. This project is designed to be runnable by anyone without commercial software licenses.

**Why representative days and not a full year?**
The rolling horizon structure means each day involves hundreds of optimization solves. Four representative days (high VER, peak demand, high forecast error, normal) demonstrate the result clearly while remaining computationally manageable on a laptop.

---

## Notes


The mathematical formulation is taken directly from the paper's Section II. No proprietary data, commercial solvers, or supercomputing resources are required. All inputs come from the publicly available RTS-GMLC test system.
