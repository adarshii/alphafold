# AlphaFold-Guided Virtual Screening with Machine Learning Rescoring

## Abstract
Structure-based virtual screening (SBVS) remains central to early drug discovery, but practical deployment is often constrained by limited access to high-resolution protein structures, uncertainty in docking-based affinity ranking, and weak reproducibility across computational environments.

This repository presents a modular, configuration-driven SBVS framework that combines protein structure prediction (AlphaFold2/ColabFold/ESMFold-compatible logic), binding pocket prioritization, ligand preparation, molecular docking, machine-learning rescoring, ADMET triage, and interactive visualization into one transparent workflow. The novelty is systems-level: each stage is isolated, auditable, and replaceable while preserving end-to-end traceability from target definition to ranked hit export.

The implementation in this repository is intentionally scaffolded for reproducible development, testing, and method communication, while the surrounding methodology supports extension to production-scale campaigns. We provide a benchmark-oriented protocol centered on SARS-CoV-2 Mpro (PDB: 6LU7), discuss expected metrics (ROC-AUC, EF<sub>1%</sub>, RMSD), and explicitly analyze known error sources including structure uncertainty, docking bias, overfitting, data leakage, and affinity-proxy mismatch.

**Keywords:** Structure-based drug discovery, molecular docking, AlphaFold, machine learning rescoring, SBVS, protein structure prediction, cheminformatics

---

## Scope and Implementation Status (Important)

This repository currently provides a **deterministic scaffold and educational reference implementation** (including dry-run mode), not a fully productionized HTVS platform. The scientific framework documented below is intentionally broader than the present mock-data execution path and is designed to support rigorous extension into real-world campaigns.

---

## 1) Scientific Background and Context

### 1.1 Structure-based virtual screening (SBVS)
SBVS evaluates large molecular libraries against a 3D target structure to prioritize candidates for biochemical follow-up. The process is conceptually simple but mathematically difficult because it combines:

1. **Pose prediction** (where/how a ligand binds).
2. **Affinity ranking** (which ligand is likely stronger).

In industrial and academic settings, SBVS is attractive because it can reduce wet-lab burden by enriching libraries before assays.

### 1.2 Why traditional SBVS pipelines are limited
Classical pipelines usually depend on experimentally solved structures (X-ray/cryo-EM), proprietary software, and heavy expert curation. Main bottlenecks include:

- **Structure availability:** many targets lack crystal structures.
- **Scoring uncertainty:** docking scores alone often correlate only moderately with experimental affinity.
- **Engineering reproducibility:** ad-hoc scripts make audit and replication difficult.

### 1.3 AlphaFold’s role in democratizing structure prediction
AlphaFold2 and derivatives (ColabFold, ESMFold) lowered the barrier to structure-guided discovery by producing high-quality structural hypotheses at scale. This is transformative for novel targets and resource-constrained labs. However, predicted structures are not guaranteed holo states and should be interpreted with confidence metrics (e.g., pLDDT), local uncertainty checks, and biological context.

### 1.4 Why docking scores are insufficient
In many benchmarks, docking-score/affinity correlations are only moderate (often reported around **r ~ 0.5–0.7**). Common causes:

- **Size bias:** larger ligands receive artificially favorable scores.
- **Lipophilicity bias:** hydrophobic compounds can be over-rewarded.
- **Entropy underestimation:** rigid receptor assumptions underrepresent conformational entropy and solvent effects.

This motivates machine-learning rescoring to learn data-driven corrections.

---

## 2) End-to-End Pipeline Architecture

### 2.1 Conceptual workflow (7 stages)

```text
Input target + YAML config
          |
          v
[1] Structure preparation/prediction
          |
          v
[2] Pocket detection + druggability scoring
          |
          v
[3] Ligand library preparation
          |
          v
[4] Docking (pose + raw score)
          |
          v
[5] ML rescoring (feature-based ranking)
          |
          v
[6] ADMET filtering
          |
          v
[7] Rank/export + dashboard visualization
```

### 2.2 Implemented stage order in code
The orchestrator (`src/alphafold_vs_pipeline/pipeline.py`) executes:

```python
structure_info = prepare_structure(config["target"])
pockets = detect_pockets(config["pocket_detection"])
compounds = prepare_library(config["library"])
docked = dock_batch(compounds, config["docking"])
rescored = rescore_poses(docked, config["rescoring"])
passed = filter_admet(rescored, config["admet"])
ranked = rank_and_export(passed, config["ranking"], out_dir)
```

This explicit order provides transparent information flow and stage-level testability.

### 2.3 Configuration-driven parameterization
All major decisions are controlled in YAML (`configs/pipeline.yaml`), including predictor, pocket tools, docking engine, and ADMET thresholds.

```yaml
target:
  name: SARS-CoV-2 Mpro
  pdb_id: 6LU7
  predictor: colabfold

pocket_detection:
  tools: [fpocket, dogsitescorer]

docking:
  engine: vina
  exhaustiveness: 8
  num_modes: 10

rescoring:
  model_type: lightgbm
```

### 2.4 Modular design principles
- Stage isolation (single responsibility)
- Composable interfaces
- Deterministic dry-run pathway for CI/testing
- Auditable outputs (`summary.json`, `hits.csv`, `hits.sdf`)

---

## 3) Quick Start (existing usage preserved)

### 3.1 Install

```bash
cd <repository_root>
pip install -e .
```

For tests:

```bash
pip install -e '.[dev]'
pytest -q
```

For dashboard support:

```bash
pip install -e '.[dashboard]'
```

### 3.2 Run demo pipeline

```bash
alphafold-vs run --config configs/pipeline.yaml --output outputs/demo --dry-run
```

### 3.3 Launch dashboard

```bash
streamlit run dashboard/app.py
```

---

## 4) Module-Wise Deep Explanation

## 4.1 Stage 1 — Protein Structure Prediction

### Scientific purpose
Rapid target characterization when experimental structures are unavailable or incomplete.

### Algorithms and model families
- **AlphaFold2:** high-accuracy deep learning structure prediction using MSAs and attention-based architecture.
- **ColabFold:** practical AF2 deployment with MMseqs2 acceleration and lower setup burden.
- **ESMFold:** protein language-model-based structure inference for fast throughput.

### Confidence calibration (pLDDT)
- **>90:** highly reliable local geometry
- **70–90:** generally reliable
- **50–70:** caution (flexible/loop regions)
- **<50:** avoid high-confidence mechanistic interpretation

### Assumptions and limitations
- Predicted states are often apo-like, not ligand-induced holo conformations.
- Side-chain and microstate uncertainty propagates into docking.
- Multi-state dynamics are incompletely represented.

### Integration in this repository
`prepare_structure(...)` records source, predictor, and structure path metadata for downstream traceability.

### Alternatives
Homology modeling, Rosetta refinement, cryo-EM model optimization, or experimental structure-first workflows.

---

## 4.2 Stage 2 — Binding Pocket Detection

### Scientific purpose
Identify and prioritize ligandable cavities to constrain docking search space.

### Geometric/chemical principles
Pocket calling integrates cavity geometry (depth, enclosure, volume) with local chemistry (hydrophobicity, H-bond capacity, charge context).

### fpocket
fpocket uses alpha-sphere geometry (related to Voronoi-like cavity characterization) for robust geometric pocket detection.

### DoGSiteScorer
DoGSiteScorer uses descriptor-driven scoring and ML-informed druggability estimates.

### Consensus strategy
Using both tools can reduce false positives versus single-tool ranking by prioritizing overlap and agreement.

### Limitations
- Cryptic/allosteric pockets may be missed in static structures.
- Results are sensitive to protonation and side-chain orientation.

---

## 4.3 Stage 3 — Ligand Library Preparation

### Scientific purpose
Transform heterogeneous sources into chemically standardized, docking-ready compounds.

### Typical processing chain
1. Active extraction from ChEMBL.
2. Decoy selection from ZINC (DUD-E style property matching).
3. SMILES canonicalization and normalization.
4. Protonation assignment at physiological pH (~7.4).
5. 3D conformer generation and minimization (e.g., MMFF94 workflows).
6. QC: duplicates, valence sanity, property range checks.

### Assumptions and limitations
- Protonation/tautomer uncertainty can dominate downstream rank shifts.
- Single-conformer docking under-samples flexible ligands.

### Alternatives
RDKit/Open Babel/OMEGA-centric preparation pipelines, expanded tautomer/protomer ensembles.

---

## 4.4 Stage 4 — Molecular Docking

### Scientific purpose
Estimate plausible binding poses and produce first-pass affinity ranking.

### Problem formulation
Docking approximates:

\[
\hat{x} = \arg\min_x E_{\text{score}}(R, L, x)
\]

where \(R\) is receptor, \(L\) ligand, and \(x\) pose/conformation variables.

### AutoDock Vina components
Vina combines stochastic/global search and local optimization with an empirical scoring approximation (sterics, hydrophobics, H-bond-related terms, conformational penalties).

### Alternative engine: Gnina
Gnina augments docking with deep-learning scoring and can improve enrichment on some targets.

### Critical parameters
- `exhaustiveness`: search depth
- `num_modes`: number of returned poses
- Grid definition: often the largest practical source of variance

### Complexity and parallelization
With fixed settings, cost scales roughly linearly with ligand count and is embarrassingly parallel across compounds.

### Post-docking QC
For benchmarked redocking, **RMSD < 2 Å** is standard near-native criterion.

---

## 4.5 Stage 5 — Machine Learning-Driven Rescoring

### Scientific purpose
Improve rank quality by learning nonlinear corrections to docking-score biases.

### Why docking fails alone
- Size and lipophilicity biases
- Weak treatment of entropy/solvation
- Force-field simplification and rigid receptor assumptions

### RF-Score feature engineering
RF-Score-style representations use protein-ligand atom-pair contact counts (commonly ~250 features depending on implementation details).

### Descriptor augmentation
Useful additions include TPSA, LogP, rotatable bonds, and pocket druggability descriptors.

### LightGBM architecture rationale
- Leaf-wise tree growth
- Strong performance on tabular descriptors
- Regularization support
- Fast inference and straightforward feature-importance analysis

### Training strategy
Typical setup: train on PDBbind with experimental affinity labels, validate on target-disjoint splits.

### Overfitting mitigation
- L1/L2 regularization
- early stopping
- cross-validation
- feature pruning

### Data-leakage prevention
- Protein-level split (avoid homolog contamination)
- Ligand scaffold split (avoid near-duplicate analog leakage)
- Time-aware split (avoid future-data contamination)

### Interpretability
SHAP analysis supports global and per-compound attribution, enabling chemically meaningful error analysis.

---

## 4.6 Stage 6 — ADMET Filtering

### Scientific purpose
Remove compounds unlikely to survive medicinal chemistry and translational filters.

### Rules and filters
- **Lipinski Rule of Five**
- **Veber rule**
- **Pfizer-style lipophilicity risk heuristics**
- **PAINS** substructure alerts

### External predictor ecosystem
- SwissADME for permeability, synthetic accessibility, medicinal chemistry flags
- pkCSM for ML-based ADMET endpoint estimates

### Typical impact
20–40% removal is common in practical virtual screening libraries (dataset-dependent).

---

## 4.7 Stage 7 — Visualization and Export

### Scientific purpose
Provide interpretable, reproducible outputs suitable for decision meetings and downstream assay planning.

### Streamlit dashboard in this repository
`dashboard/app.py` shows:
- run metadata summary,
- 3D scatter (docking score vs ML score),
- downloadable JSON hit export.

### Output formats
- CSV (`hits.csv`)
- SDF (`hits.sdf`)
- JSON summary (`summary.json`)

Extended workflows may add ROC plots, calibration plots, and assay-prioritization reports.

---

## 5) Machine Learning Deep Dive (Comprehensive)

### 5.1 Why docking score correlation is weak
Across many docking benchmarks, score-affinity correlation is often moderate (commonly near **r ~ 0.5–0.6** for difficult sets). This is insufficient for reliable top-fraction triage in noisy real-world campaigns.

### 5.2 Systematic scoring biases
1. **Size bias:** larger ligands tend to gain favorable contact terms.
2. **Lipophilicity bias:** hydrophobic packing can dominate score terms.
3. **Entropy omission:** conformational/solvent penalties are underrepresented.

### 5.3 ML strategy
Learn a mapping \(g(\mathbf{f}) \to \hat{y}\) from engineered descriptors \(\mathbf{f}\) to affinity proxy or active-likelihood signal, rather than relying on raw docking score alone.

### 5.4 Feature engineering details
- RF-Score contact vectors (atom-type interaction counts)
- Ligand descriptors (TPSA, cLogP, MW, rotatable bonds)
- Pocket descriptors (volume, enclosure, druggability)
- Docking descriptors (best score, pose spread, mode stability)

### 5.5 Why LightGBM
- Handles nonlinear interactions in mixed descriptor spaces
- Efficient on CPU and scalable to large tabular datasets
- Interpretability via feature importance + SHAP

### 5.6 Critical ML failure modes

#### Overfitting
Warning pattern:
- train \(R^2 = 0.85–0.90\)
- validation \(R^2 = 0.55–0.65\)
- gap \(>0.15\) is often problematic

Mitigation: stronger regularization, depth control, feature selection, repeated CV.

#### Data leakage
- Ligand-level analog leakage
- Protein-level homology leakage
- Temporal leakage (future complexes in training)

Mitigation: scaffold/protein/time-aware splits aligned to deployment.

#### Domain shift
Training data (e.g., PDBbind) is biased toward specific target classes (kinases/proteases). Out-of-domain targets (e.g., some GPCRs/membrane proteins) can degrade sharply.

#### Affinity proxy mismatch
IC<sub>50</sub> is not equivalent to \(\Delta G_{bind}\); assay kinetics, transport, and cellular context confound equilibrium affinity interpretations.

### 5.7 Interpretability and decision support
SHAP values can explain individual rank changes, helping prioritize compounds with mechanistically plausible signals over artifacts.

---

## 6) Benchmark Dataset: SARS-CoV-2 Mpro (PDB 6LU7)

### 6.1 Selection rationale
Mpro is a practical benchmark because it has strong structural data, extensive inhibitor literature, and broad community familiarity.

### 6.2 Recommended benchmark composition
- **Actives:** 157 compounds (ChEMBL, IC<sub>50</sub> < 10 µM)
- **Decoys:** 10,000 ZINC15 compounds (DUD-E-like matching: MW ±50 Da, LogP ±0.5)
- **Class balance:** ~1.5% actives vs 98.5% decoys (realistic screening imbalance)

### 6.3 Biases and limitations
- Publication bias toward positives
- Protease-scaffold concentration
- Decoy-construction artifacts in DUD-E-style sets

### 6.4 Mitigation
External validation with independent sources (BindingDB/literature curation), scaffold-aware analysis, and prospective assay confirmation.

---

## 7) Evaluation and Validation Strategy

### 7.1 ROC-AUC
Measures threshold-independent ranking discrimination:
- 0.5 random
- 0.7+ useful
- 0.8+ strong (context dependent)

### 7.2 Enrichment Factor (EF<sub>1%</sub>)
\[
EF_{x\%} = \frac{\text{actives in top }x\% / N_{x\%}}{\text{total actives} / N_{\text{all}}}
\]

Interpretation:
- EF<sub>1%</sub> > 5: good
- EF<sub>1%</sub> > 10: excellent (dataset/protocol dependent)

### 7.3 RMSD-based pose accuracy
Near-native threshold: **RMSD < 2 Å**. Typical successful docking campaigns can reach ~70–85% redocking success under appropriate setup.

### 7.4 Calibration curves
Compare predicted and observed activity frequencies across score bins; essential when translating model outputs into selection quotas.

### 7.5 Train-test design
- 80/20 split by protein/scaffold (not naive random ligand split)
- 5-fold stratified CV
- report mean ± standard deviation

### 7.6 Prospective validation protocol
1. Train on historical complexes (e.g., PDBbind).
2. Predict on Mpro virtual library.
3. Experimentally test top-ranked compounds.
4. Report ROC-AUC, EF<sub>1%</sub>, and hit rate.

### 7.7 Expected outcomes on Mpro benchmark (illustrative)
- Vina ROC-AUC: **0.72 ± 0.05**
- LightGBM ROC-AUC: **0.80 ± 0.05** (Δ +0.08)
- Vina EF<sub>1%</sub>: **4.2×**
- LightGBM EF<sub>1%</sub>: **6.5×** (≈ +55%)

---

## 8) Critical Limitations (Honest Assessment)

1. **Protein structure assumptions**
   - Static models cannot represent full conformational ensembles.
   - Low-confidence regions (pLDDT < 50) are unreliable for detailed docking interpretation.

2. **Docking inaccuracies**
   - Decoys can outscore true binders.
   - Structured water and receptor flexibility are incompletely modeled.

3. **No direct experimental validation in this repository scaffold**
   - In silico ranking is triage, not confirmation.
   - Biochemical/cellular follow-up is mandatory.

4. **ML generalization risk**
   - Target-class shift can reduce performance sharply.

5. **Affinity proxy gap**
   - IC<sub>50</sub> ≠ \(\Delta G_{bind}\) under many assay conditions.

6. **Orthosteric bias**
   - Default setup may miss allosteric opportunities.

7. **Compute/resource barriers**
   - 10k-scale campaigns may require ~10–18 GPU-hours depending on protocol.

---

## 9) Future Work and Extensions

### 9.1 Molecular dynamics integration
Ensemble docking over MD snapshots (e.g., 100 ns, sampled every 5 ns) to account for receptor plasticity.

### 9.2 Free-energy refinement
- **MM-PBSA:** moderate accuracy, fast triage-level re-ranking.
- **FEP:** high-accuracy potential at substantially higher computational cost.

### 9.3 Multi-omics integration
Incorporate transcriptomics/proteomics context for target engagement and pathway-aware prioritization.

### 9.4 Active learning loops
Iterative predict → experimentally test uncertainty-prioritized compounds → retrain.

### 9.5 Allosteric/cryptic pocket modeling
Add dynamic pocket discovery and non-orthosteric intervention workflows.

---

## 10) Reproducibility and Engineering Design

### 10.1 Stage modularity in repository
- `src/alphafold_vs_pipeline/stages/structure.py`
- `src/alphafold_vs_pipeline/stages/pockets.py`
- `src/alphafold_vs_pipeline/stages/libraries.py`
- `src/alphafold_vs_pipeline/stages/docking.py`
- `src/alphafold_vs_pipeline/stages/rescoring.py`
- `src/alphafold_vs_pipeline/stages/admet.py`
- `src/alphafold_vs_pipeline/stages/ranking.py`

### 10.2 Config as experimental contract
All key hyperparameters are in YAML and version-controlled for auditability.

### 10.3 Deterministic execution
Production deployments should enforce fixed seeds (`random.seed(42)`, `np.random.seed(42)`, `torch.manual_seed(42)`) and log software/hardware metadata.

### 10.4 Scalability profile
- Single-machine: ~10k compounds in 10–18 GPU-hours (target/protocol dependent).
- Cloud: batch scaling for 100k compounds.

### 10.5 Dependency management
Containerized environments (Docker/Apptainer) are recommended for cross-system reproducibility.

---

## 11) Comparative Analysis: Traditional vs ML-Augmented SBVS

| Aspect | Traditional SBVS | This Pipeline (ML-Augmented) |
|---|---|---|
| Structure source | Experimental (X-ray, cryo-EM) | Predicted (AlphaFold2/ColabFold/ESMFold) + optional experimental |
| Docking engine | AutoDock/FlexX/Glide variants | AutoDock Vina (Gnina-compatible extension path) |
| Affinity prediction | Docking score only | Docking + ML rescoring (LightGBM) |
| Enrichment (EF<sub>1%</sub>) | ~3–5× | ~5–10× target range |
| Turnaround time | Weeks–months (if structure unavailable) | Hours to ~1 day for triage |
| Cost profile | High (infrastructure/licensing/structure campaigns) | Lower barrier, open-source centric |
| Accessibility | Mostly large pharma/core facilities | Academic labs/startups/biotech |
| Interpretability | High physics intuition | Hybrid (physics + SHAP-based ML explanations) |
| Generalization | Often robust with quality structures | Target-dependent; retraining advised |

**Use traditional pipelines when:** high-budget, time-critical, highly regulated programs demand maximal structure certainty.  
**Use this pipeline when:** rapid iteration, novel targets, and budget-constrained exploratory campaigns are priorities.

---

## 12) Expected Results and Typical Outputs

### 12.1 Illustrative Mpro benchmark expectations

**Docking stage**
- RMSD < 2 Å success for known binders: 70–95% (protocol dependent)
- Mean Vina score (actives): ~−8.2 ± 0.8 kcal/mol
- Mean Vina score (decoys): ~−6.5 ± 1.2 kcal/mol

**Rescoring stage**
- ROC-AUC (Vina): 0.72 ± 0.05
- ROC-AUC (LightGBM): 0.80 ± 0.05
- EF<sub>1%</sub> (Vina): 4.2×
- EF<sub>1%</sub> (LightGBM): 6.5×

**ADMET filtering**
- Expected removal: 20–40%
- Residual library from 10,000 compounds: ~6,000–8,000

### 12.2 Output artifacts

```text
outputs/demo/
├── summary.json
├── hits.csv
├── hits.sdf
└── logs/           # optional in extended workflows
```

---

## 13) Practical Guidance and Troubleshooting

### 13.1 Missing Plotly/Streamlit modules
Install dashboard extras:

```bash
pip install -e '.[dashboard]'
```

### 13.2 Dashboard says summary file is missing
Run:

```bash
alphafold-vs run --config configs/pipeline.yaml --output outputs/demo --dry-run
streamlit run dashboard/app.py
```

### 13.3 Unstable rankings across runs
- Fix seeds
- Pin dependencies
- Keep config and structure inputs identical
- Log environment metadata

### 13.4 Poor enrichment
- Re-check pocket box definition
- Audit protonation/tautomer generation
- Test leakage-safe split strategy
- Recalibrate/retrain model per target class

---

## 14) Figure Suggestions for Publication/Portfolio

1. **Pipeline architecture** (7-stage flow + config layer)
2. **ROC comparison** (Vina vs ML rescoring, mean ± CI)
3. **Early enrichment plot** (EF at 0.5/1/2/5%)
4. **Calibration curve** (predicted vs observed activity)
5. **SHAP summary** (global feature attribution)

### ASCII ROC sketch

```text
TPR
1.0 |                             * ML-rescored
    |                         *
0.8 |                     *
    |                 *
0.6 |             *          baseline docking
    |         *
0.4 |      *
    |   *
0.2 | *
    +------------------------------------------- FPR
      0.0       0.2       0.4       0.6      1.0
```

---

## 15) Insight Notes

> **Insight 1:** Docking is best treated as a noisy ranking prior, not a final affinity estimator.

> **Insight 2:** In ML-based rescoring, split design (leakage control) is often more important than model brand choice.

> **Insight 3:** For medicinal chemistry triage, early enrichment and calibration are usually more actionable than global regression fit.

---

## 16) Recommended Reading and References

### SBVS and benchmarking
- Leung SC, Bodkin M, von Delft F, Fink A, Morris GM. *J Chem Inf Model.* 2021;61(5):1886–1896.
- Huang N, Shoichet BK, Irwin JJ. *J Chem Inf Model.* 2006;46(1):243–253.
- Morris GM, Huey R, Lindstrom W, et al. *J Comput Chem.* 2009;30(16):2785–2791.

### Protein structure prediction
- Jumper J, Evans R, Pritzel A, et al. *Nature.* 2021;596(7873):583–589.
- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M. *Nat Methods.* 2022;19(6):679–682.
- Lin Z, Akin H, Rao R, et al. *bioRxiv.* 2023.

### Machine learning in drug discovery
- Ballester PJ, Mitchell-Koch J. *J Chem Inf Model.* 2010;50(5):816–834.
- Ke G, Meng Q, Finley T, et al. *NeurIPS.* 2017.
- Ragoza M, Hochuli J, Idrobo E, Sunseri J, Koes DR. *J Chem Inf Model.* 2017;57(4):942–957.

### ADMET and drug-likeness
- Daina A, Michielin O, Zoete V. *Sci Rep.* 2015;5:13230.
- Pires DEV, Blundell TL, Ascher DB. *J Chem Inf Model.* 2015;55(10):1737–1747.
- Baell JB, Holloway GA. *J Med Chem.* 2010;53(7):2719–2740.

### SARS-CoV-2 context
- Liu X, Wang XJ. *Viruses.* 2020;12(4):406.

---

## 17) Extended Methodological Notes

### 17.1 Stage-by-stage I/O contract (practical engineering view)

| Stage | Input | Core transformation | Output | Failure modes to monitor |
|---|---|---|---|---|
| Structure | Target metadata, sequence/structure source | Prediction or structure registration | Receptor structure + confidence metadata | Missing residues, low-confidence loops |
| Pocketing | Receptor coordinates | Cavity detection + druggability scoring | Ranked pocket list + geometric descriptors | False-positive cavities, missed cryptic pockets |
| Library prep | SMILES sets (actives/decoys) | Standardization, protonation, 3D conformers | Docking-ready ligand set | Invalid valence, tautomer drift, duplicates |
| Docking | Receptor + ligands + search space | Pose search + empirical scoring | Pose ensemble + raw docking scores | Grid misplacement, unrealistic poses |
| ML rescoring | Docked complexes + descriptors | Feature extraction + model inference | Re-ranked compounds | Leakage, overfitting, domain shift |
| ADMET | Ranked compounds | Rule-based + predictive filtering | Developability-filtered shortlist | False negatives on novel chemotypes |
| Export/UI | Filtered hits + metadata | Serialization + interactive visualization | CSV/SDF/JSON + dashboard assets | Incomplete metadata, dashboard dependency issues |

This table is useful when converting this scaffold into production-grade screening infrastructure: each stage can be validated independently, then integrated with explicit contracts.

### 17.2 Docking algorithmic detail and search assumptions
Although practical pipelines often describe docking as a single black-box command, two independent components govern quality:

1. **Search procedure** (global + local optimization over pose/conformation space).
2. **Scoring approximation** (energy-like surrogate for ranking).

AutoDock-family methods historically combine stochastic population-based exploration and local refinement (often discussed in the context of Lamarckian ideas in AutoDock4 literature), while Vina emphasizes efficient iterative optimization under an empirical score. In both cases, exhaustive physics is not solved exactly; practical docking is a heuristic optimization under strict compute constraints. Therefore:

- better search does not guarantee better ranking if the scoring surrogate is biased,
- better scoring terms can still fail if the search misses near-native modes.

For this reason, quality control should include both **pose plausibility checks** and **rank-performance checks**.

### 17.3 On score decomposition and why ranking errors persist
A simplified empirical score can be conceptualized as:

\[
S \approx w_1 E_{\text{vdW}} + w_2 E_{\text{electrostatic}} + w_3 E_{\text{H-bond}} + w_4 E_{\text{hydrophobic}} + w_5 E_{\text{torsion penalty}}
\]

Even if this decomposition captures coarse trends, persistent mismatch arises because:

- receptor flexibility is often frozen,
- structured waters are omitted or approximated,
- long-range polarization and context-dependent protonation are simplified,
- entropy and kinetics are incompletely represented.

Hence, top-ranked docking compounds are best interpreted as **enriched hypotheses**, not direct surrogates for biochemical truth.

### 17.4 ML model governance checklist (for publication-grade rigor)
When converting rescoring experiments into claims suitable for methods papers:

1. Define the prediction target explicitly (classification vs regression, endpoint definitions).
2. Declare split strategy before model training.
3. Report confidence intervals, not only point metrics.
4. Audit feature leakage (especially target-derived or assay-derived proxies).
5. Provide negative controls (e.g., shuffled labels, random baselines).
6. Publish failure cases (targets where rescoring underperformed docking).

This governance framework is often more valuable than small metric gains because it determines scientific credibility.

### 17.5 Data leakage in practice: concrete examples
- **Ligand leakage:** same Bemis–Murcko scaffold appears in both train and test; model memorizes chemotype motifs.
- **Protein leakage:** homologous binding pockets share interaction signatures across split folds.
- **Temporal leakage:** complexes deposited after a given date appear in train while benchmarking “future” targets.

In all three scenarios, apparent ROC-AUC improvements can be inflated. Recommended mitigation includes scaffold-aware splitting, sequence-similarity thresholding for proteins, and strict chronological holdouts.

### 17.6 Domain shift diagnostics
A model performing well on kinases/proteases can fail on GPCRs/membrane proteins because pocket environments, hydration profiles, and ligand chemotypes differ. Practical diagnostics:

- monitor calibration drift by target family,
- compare feature distributions between train and inference sets,
- compute per-family enrichment metrics rather than only global averages.

If distribution mismatch is large, perform transfer calibration or retraining.

### 17.7 Affinity proxy caveat (IC50, Ki, and \(\Delta G\))
For an educationally honest pipeline narrative, distinguish:

- **\(K_i\)**: thermodynamic inhibition constant (closer to equilibrium affinity in idealized settings),
- **IC\(_{50}\)**: assay-condition-dependent concentration reducing activity by 50%,
- **\(\Delta G_{bind}\)**: free energy change under defined thermodynamic model.

These are related but not interchangeable without assumptions. Reporting this distinction prevents over-interpretation of ML “affinity” outputs.

### 17.8 Recommended reporting template for results sections
For each benchmark target, report:

1. Dataset curation criteria and exclusion counts.
2. Pocket protocol and docking grid specification.
3. Docking baseline metrics (AUC, EF<sub>1%</sub>, pose RMSD).
4. Rescoring metrics with uncertainty.
5. ADMET attrition statistics.
6. Top-k hit composition by scaffold family.
7. Error analysis of false positives/false negatives.

This structure aligns with review expectations in computational medicinal chemistry journals.

---

## 18) Benchmark-Specific Practical Expectations (SARS-CoV-2 Mpro)

### 18.1 Target rationale and structural context
SARS-CoV-2 Mpro is suitable for pedagogical and methodological benchmarking because:

- multiple high-quality structures are available,
- inhibitor chemistry is rich and diverse enough for statistical analysis,
- the target sits at the intersection of virology and medicinal chemistry, enabling cross-domain communication.

The representative structure **6LU7** (resolution 1.88 Å) is a widely used reference for docking methodology demonstrations.

### 18.2 Suggested dataset composition

| Component | Value | Notes |
|---|---|---|
| Actives | 157 | ChEMBL-derived, IC<sub>50</sub> < 10 µM criterion |
| Decoys | 10,000 | ZINC15, DUD-E style property matching |
| Class balance | 1.5% active / 98.5% decoy | Reflects realistic screening imbalance |

### 18.3 Decoy matching assumptions
Common DUD-E-style matching constraints include property windows such as MW ±50 Da and LogP ±0.5. These constraints improve fairness versus random decoys but can still introduce synthetic benchmark artifacts.

### 18.4 Expected stage-wise outcomes (illustrative, not guaranteed)

| Stage | Metric | Typical target value |
|---|---|---|
| Docking baseline | ROC-AUC | 0.72 ± 0.05 |
| ML rescoring | ROC-AUC | 0.80 ± 0.05 |
| Docking baseline | EF<sub>1%</sub> | 4.2× |
| ML rescoring | EF<sub>1%</sub> | 6.5× |
| Pose quality | Redocking RMSD < 2 Å | 70–85% typical range |
| ADMET filtering | Removal rate | ~20–40% |

### 18.5 Top-fraction interpretation example
In a 10,000-compound campaign with 1.5% actives, random selection of top 1% (100 compounds) would expect ~1.5 actives. EF<sub>1%</sub> = 6.5 implies expected active recovery near 9–10 compounds in top 100 (subject to dataset and assay variation).

### 18.6 Prospective validation protocol (recommended)
1. Train rescoring on historical data only.
2. Freeze model and preprocessing.
3. Run blind prediction on benchmark library.
4. Select top 100 for biochemical assay.
5. Compute experimental hit rate and compare against docking-only baseline.
6. If possible, expand top 10 into orthogonal cellular and selectivity assays.

This protocol aligns computational metrics with translational decision points.

---

## 19) Reproducibility and Deployment Blueprint

### 19.1 Minimal reproducibility checklist
- Pin dependency versions (Python + core libraries).
- Version control config files and datasets.
- Store command lines and random seeds.
- Archive intermediate artifacts (prepared structures, pocket metadata, docking poses, feature matrices).
- Record hardware profile (CPU/GPU type and driver stack).

### 19.2 Deterministic execution guidance
In production codepaths, explicitly set deterministic controls:

```python
import random
import numpy as np

random.seed(42)
np.random.seed(42)
# torch.manual_seed(42) when torch is in the runtime stack
```

Deterministic settings do not eliminate all floating-point variability but significantly reduce irreproducible rank drift.

### 19.3 Container and environment strategy
Use Docker/Apptainer images to lock:

- OS and compiler stack
- BLAS/LAPACK backends
- Python package versions
- GPU runtime (if applicable)

This avoids hidden environment-dependent variation across labs and cloud nodes.

### 19.4 Scalability profile
- **Single machine:** ~10k compounds in ~10–18 GPU-hours (target and settings dependent).
- **Cloud mode:** batch partitioning enables 100k+ compound campaigns by horizontal scaling.

For cloud reproducibility, preserve immutable image tags and parameter snapshots.

### 19.5 Suggested repository extension map
The current scaffold can be extended by replacing stage internals while keeping interfaces stable:

- `structure.py`: integrate actual AlphaFold/ESMFold runners and confidence parsers.
- `pockets.py`: plug fpocket/DoGSiteScorer wrappers with standardized schema.
- `libraries.py`: integrate RDKit/Open Babel canonicalization and conformer generation.
- `docking.py`: connect Vina/Gnina execution engine.
- `rescoring.py`: attach trained LightGBM models + feature extraction.
- `admet.py`: include full ADMET predictors and PAINS SMARTS filtering.
- `ranking.py`: extend report generation and publication-ready figures.

This interface-preserving strategy supports iterative development without breaking orchestration logic.

### 19.6 Practical publication note
For preprint or thesis inclusion, clearly distinguish:

1. **Current repository implementation status** (scaffolded deterministic flow),
2. **Methodological blueprint** (full scientific protocol),
3. **Validated experimental evidence** (if/when assays are completed).

This separation improves transparency and strengthens scientific credibility.

---

## 20) Repository Citation

Please use `CITATION.cff` for software citation metadata.

---

## 21) Scope Note

See **Scope and Implementation Status (Important)** near the top of this document.
