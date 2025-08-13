# Lipid Rendering Docking Pipeline Development Roadmap

Use this file (`ROADMAP.md`) as the definitive phase-by-phase plan. Cursor should complete each phase fully—including all tests—then mark tasks as done in `docs/current_feature.md` before moving on.

---

## Phase 0 — Enforce REAL_DOCKING & Deprecate Mock by Default ✅ COMPLETED

**Goal**: Always use AutoDock Vina unless explicitly overridden. Mock docking only if `DOCKING_ALLOW_MOCK=true`.

- [x] Add `DOCKING_ALLOW_MOCK` config (default: `false`)
  - If false and Vina unavailable → API returns HTTP 503 with clear message.
  - Implement `/api/dock/capabilities` to return `{ vina_available, vina_version, engine_default }`.
- [x] In results JSON and DB job record, always include:
  - `engine`: `"vina" | "mock"`
  - `is_mock`: boolean
- [ ] **Frontend**:
  - Display engine badge: `"AutoDock Vina vX.Y.Z"` when real,
    or red `"MOCK (demo only)"` badge when mock used.
  - If mock, disable or gray out "export pose" buttons.
- [x] **Tests**:
  - Unit test for `capabilities` endpoint with/without Vina.
  - Integration test: when `DOCKING_ALLOW_MOCK=false` & Vina missing → `/dock/run` yields 503.
  - Frontend UI snapshot: engine badge correctness.

---

## Phase 1 — Seeded & Deterministic Vina ✅

**Goal**: Same inputs + same seed yield same outputs; reproducibility guaranteed.

- [x] API: accept optional integer `seed` parameter.
  - Default: either auto-random (nondeterministic) or fixed (e.g., `42`).
  - Use in Python API: `Vina(seed=seed)` :contentReference[oaicite:1]{index=1}.
- [x] Reflect `seed` in result parameters.
- [x] Add structured logging:
  - On start/end: log `job_id`, ligand/receptor IDs, center/size, exhaustiveness, num_modes, seed, time elapsed, best affinity.
- [x] **Tests**:
  - **Deterministic**: same inputs & seed → identical energies & pose coordinates.
  - **Nondeterministic**: different seeds → at least one difference.
  - Parameter validation: appropriate ranges and integer-only seed.

---

## Phase 2 — True Docked Poses to Frontend ✅

**Goal**: Scene rendered using real pose geometries, not translations.

### Backend
- [x] After `v.dock(...)`, extract poses using:
  - `v.write_poses("out.pdbqt", n_poses=...)`, or
  - `v.poses()` :contentReference[oaicite:2]{index=2}.
- [x] Convert each pose to SDF using **Meeko** (preferred for fidelity) or OpenBabel as fallback with warning :contentReference[oaicite:3]{index=3}.
- [x] In JSON results, each pose must include `{ mode, affinity, rmsd_lb, rmsd_ub, sdf }`.

### Frontend
- [ ] Replace translation hack: load each `pose.sdf` directly:
  - Use `viewer.addModel(sdf, 'sdf')` (3Dmol.js supports SDF) :contentReference[oaicite:4]{index=4}.
- [ ] Validate spatial differences between poses; use distinct colors/styles.

### Tests
- [x] Backend: SDF exists per pose; valid RDKit parse.
- [x] Coordinates differ across modes.
- [ ] Frontend: snapshot tests or simple e2e checks confirm positional differences on UI.

---

## Phase 3 — Deterministic Preparation (Ligand & Receptor) ✅

**Goal**: Input preprocessing fully deterministic and documented.

### Ligand
- [x] Generate 3D using RDKit ETKDG with fixed `randomSeed`.
- [x] Add hydrogens; optionally MMFF minimize.
- [x] Document any protonation/tautomer policy.

### Receptor
- [x] Clean PDB, add hydrogens, assign charges.
- [x] Use **Meeko** for PDBQT prep :contentReference[oaicite:5]{index=5}.

### Tests
- [x] SMILES → 3D → PDBQT → dock → SDF pipeline yields repeatable output.
- [x] Snapshot tests: atom count, charge consistency, etc.

---

## Phase 4 — UI Visibility & Metadata ✅

**Goal**: Transparency in parameters and engine details for the user.

- [x] Show in job detail panel:
  - Engine, Vina version, seed, time runtime, exhaustiveness, num_modes, center/size.
- [x] Result header badges: green for real, red for mock.
- [x] Export includes metadata JSON alongside pose SDF.

### Tests
- [x] API returns all metadata correctly.
- [x] Frontend rendering of metadata and badges (snapshot).
- [x] Accessibility: ARIA labels & color-blind friendly indicators.

---

## Phase 5 — CI with Seeded Smoke Test ✅

**Goal**: Every PR builds, tests, runs a fast seeded docking.

- [x] GitHub Actions workflow:
  - Setup micromamba/conda (cached) with Vina, RDKit, Meeko.
  - Run unit tests (backend + frontend).
  - Run a fast seeded docking (small sample, e.g. exhaustiveness=4, modes=3, seed=42).
  - Persist logs, JSON, and SDF as artifacts.
- [x] (Optional) Docker build step to validate production container.

### Tests
- [x] CI passes clean build with all tests.
- [x] Seeded smoke test outputs identical across runs.

---

## Phase 6 — Benchmark Suite & Regression Tests ✅

**Goal**: Empirical validation; detect regression via known docking cases.

**Select benchmarks** (use known PDB complexes):
- Streptavidin–Biotin (`1STP`)
- HIV-1 Protease + Inhibitor (`1HVR`)
- DHFR + Methotrexate (`4DFR`)
- Trypsin + Benzamidine (`3PTB`)
- Estrogen Receptor + Estradiol (`1ERE`)

**Tasks**
- [x] Add `benchmarks/` folder: assets or scripts to fetch & prep PDBs.
- [x] Create seed-based docking scripts: compute RMSD vs. crystal (via RDKit align), record affinity.
- [x] Set thresholds: e.g., top-pose RMSD ≤ 2 Å; energy variance ≤ 0.2 kcal/mol.
- [x] CI:
  - Run lightweight subset on all PRs.
  - Full set nightly or on release tags.

---

## Phase 7 — Documentation (Developer & User Guides) ✅

- [x] `DOCKING_MODES.md`: real vs. mock policy; default behavior; mock indicator.
- [x] `REPRODUCIBILITY.md`: explain seeds, deterministic embeddings, parameter logging.
- [x] `VISUALIZATION.md`: pose extraction (Vina APIs), Meeko usage, 3Dmol rendering.
- [x] `CI.md`: explain GitHub Actions, fast smoke, nightly benchmarks.
- [x] Update `README.md`: badges for build status and "Uses AutoDock Vina" messaging.
- [x] `CONTRIBUTING.md`: include testing rules (tests required per feature; smoke must pass).

---

## Optional Phase 8 — Performance Tuning & Polishing

- [ ] Expose `cpu` parameter for Vina to tune threading :contentReference[oaicite:7]{index=7}.
- [ ] Structured logging (JSON logs) and internal `/metrics` dashboard.
- [ ] Preset docking configs (“fast”, “balanced”, “thorough”) mapped to parameter sets.

---

## Definition of Done

- [ ] Vina is the **default docking engine**; mock only if explicitly enabled.
- [ ] Seeded runs reproducible; documented and tested.
- [ ] True docked poses visually represented; no translation hacks.
- [ ] All CI workflows operational: unit tests, smoke test, optional benchmarks, Docker build.
- [ ] Benchmark suite provides reproducible performance validation.
- [ ] Comprehensive documentation for developers and users.

---

### Testing References

- **Vina seed parameter**: `Vina(seed=...)`, check via `Vina(...)` class docs :contentReference[oaicite:8]{index=8}.
- **Reproducibility note**: manual enforces same seed yields exact results :contentReference[oaicite:9]{index=9}.
- **3Dmol viewer supports SDF**: via `addModel(data, "sdf")` :contentReference[oaicite:10]{index=10}.

---

Once each phase is complete, **remember**:
- Update `docs/current_feature.md` by marking the relevant tasks as `[x]` before proceeding to the next phase.
- Run **all tests** again to confirm clean state.
- Then, continue to the next phase.

Good luck, Cursor!
::contentReference[oaicite:11]{index=11}
