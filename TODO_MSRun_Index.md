TODO: Use Full-Proteome Index For Wrong IDs In MSRun

Summary
- Goal: generate wrong peptide identifications in MSRun from the full proteome, not just the observed peptide table.
- Constraint: do NOT put the index into `parameters`; keep it small. Pass the index as an explicit function argument where needed.

API Changes (minimal and explicit)
- Build the index separately and pass it to MSRun as an additional argument:
  - Change signature to: `MSRunSim(Digested, parameters, search_index = NULL, use_index_for_wrong_ids = TRUE, wrong_id_mass_tol_ppm = 20)`.
  - Keep defaults so existing calls continue to work (no breaking change).

Where To Build The Index
- Build the index over the full proteome BEFORE expression filtering. Two options:
  1) Call `build_search_index_from_df(all_proteins_df, parameters)` right after FASTA cleaning (preferred, already available).
  2) Or use the existing standalone `build_search_index(parameters)` with `PathToFasta` and matching digestion settings.
- Do NOT store the index in `parameters`. Instead:
  - Capture it in `run_sims()` after ground-truth creation and pass it as `search_index` into `MSRunSim()` calls.
  - Alternatively, return it from the ground-truth stage alongside proteoforms, but pass it explicitly, not nested into `parameters`.

Implementation Steps
1) Pipeline plumbing
   - In `run_sims()` (R/00_BatchRunFuncs.R):
     - After sample preparation (or protein import), build or retrieve the full-proteome index once per run.
     - Store it in a local variable (not in `Param`) and pass it to each `MSRunSim()` call via the new `search_index` argument.

2) MSRun wrong-ID assignment (R/03_MSRun.R)
   - Compute `n_wrong <- round(parameters$WrongIDs * nrow(MSRun))`.
   - If `use_index_for_wrong_ids && !is.null(search_index) && n_wrong > 0`:
     - Sample proteins by index weight: `pi <- sample(seq_along(search_index$proteins), n_wrong, replace = TRUE, prob = search_index$weights)`.
     - Cache peptide templates per unique `pi` using `protein_peptide_template_from_index()`.
     - For each selected MSRun row i:
       - Optionally filter the cached template by mass tolerance: `abs(1e6*(template$MZ1 - MSRun$MZ1[i])/MSRun$MZ1[i]) <= wrong_id_mass_tol_ppm`.
       - If filtered set is empty, fall back to unfiltered template; if still empty, fallback to current shuffle method.
       - Pick one decoy peptide; overwrite identification fields on row i:
         - `Accession`, `Peptide`, `Start`, `Stop`, `MC` from decoy.
         - Reset PTM annotations: `PTMPos <- list()` and `PTMType <- list()` for that row.
         - Set `WrongID <- TRUE`.
     - Keep intensities and detectability unchanged.
   - Else: fallback to current behavior (shuffle intensities among selected rows) and set `WrongID` accordingly.

3) Parameters (no bloat)
   - Keep existing `WrongIDs` and `WrongLocalizations`.
   - Add optional MSRun args only (NOT in `parameters`):
     - `search_index = NULL`, `use_index_for_wrong_ids = TRUE`, `wrong_id_mass_tol_ppm = 20`.

4) Performance considerations
   - Template caching: build per-protein peptide templates once per `MSRunSim()` call.
   - If needed later, augment the index to carry precomputed templates to remove runtime expansion.

5) Testing
   - Unit test that with a small `WrongIDs` and a provided `search_index`, output rows flagged `WrongID` have:
     - Identification fields replaced; intensities unchanged.
     - `WrongID` count equals `round(WrongIDs * nrow(MSRun))` (within bounds if corner cases).
     - When `wrong_id_mass_tol_ppm` is small, selected decoys are within tolerance, with fallback exercised when needed.
   - Backward-compat: with `search_index = NULL`, behavior matches existing shuffle logic.

6) Rollout & Cleanups
   - Wire index build in `run_sims()` after protein input to ensure full set coverage.
   - Do not store index in `parameters`; pass explicitly into `MSRunSim()` only.
   - Optionally remove returning the index from `proteinInput()` once plumbing is complete to avoid duplication.

Notes
- The index is built over the fully cleaned FASTA (unusual AAs and duplicate accessions removed), but BEFORE `PercExpressedProt` filtering, as requested.
- Digestion can remain purely `fastDigest` for now; the index is only used for wrong-ID sampling in MSRun.

