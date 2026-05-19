"""
TROTS Proton Optimizer Comparison: C++ IPOPT vs MATLAB/matRad IPOPT
--------------------------------------------------------------------
Both solvers minimize the same objective:
    f(w) = sum_i  penalty_i * sum(max(+/-(D_i*w - bound_i), 0)^2) / n_i
with hard inequality constraints via logsumexp approximation (eps=CONSTRAINT_EPS).

Usage
-----
1. For each case you want to compare against MATLAB, run the TROTS patient in
   matRad using readTROTSPatient.m (see
   https://github.com/SebastiaanBreedveld/TROTS/blob/master/matRad_import/readTROTSPatient.m).
   The script automatically saves the optimized weights to Matlab_ProtonXX_w.txt
   in the patient folder; point MATLAB_W_DIR at that folder.
   If no weight file is present for a case the MATLAB comparison is skipped.

2. Set CASE to the patient number you want to run (1–20).

3. Run inside 3D Slicer's (version > XX.YY) Python console (requires the SlicerRT
   ExternalBeamPlanning module to be loaded for the C++ IPOPT solver):
       exec(open('/path/to/compare_vs_matlab.py').read())

   All config variables (CASE, TROTS_DIR, MATLAB_W_DIR, OUTPUT_DIR, …) can be
   set as globals before calling exec() — the script picks up pre-defined values
   and falls back to the defaults in the Configuration section otherwise.
   You can also edit the defaults directly in the section below and just call exec().
       CASE = 3
       TROTS_DIR    = r'E:\data\TROTS\Protons'
       MATLAB_W_DIR = r'E:\data\TROTS\Protons'
       OUTPUT_DIR   = r'E:\data\TROTS\output'
       exec(open('/path/to/compare_vs_matlab.py').read())

Outputs (written to OUTPUT_DIR)
--------------------------------
   dvh_comparison_TROTSProtonXX.png  — per-structure DVH plots + objective table
   dvh_comparison_TROTSProtonXX.txt  — mean-dose table and constraint-violation summary
"""

import time
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sps

# ─── Configuration ────────────────────────────────────────────────────────────

CASE              = globals().get('CASE',              1)      # 1–20 → Proton 01 … Proton 20

# Optimizer settings
MAX_ITER              = globals().get('MAX_ITER',              10000)
LBFGS_HISTORY         = globals().get('LBFGS_HISTORY',         50)
ACCEPTABLE_TOL        = globals().get('ACCEPTABLE_TOL',        1e10)   # set high to disable NLP-error sub-criterion (matches matRad)
ACCEPTABLE_OBJ_CHANGE = globals().get('ACCEPTABLE_OBJ_CHANGE', 1e-4)
LINEAR_SOLVER         = globals().get('LINEAR_SOLVER',         'mumps') # 'pardisomkl'
CONSTRAINT_EPS        = globals().get('CONSTRAINT_EPS',        1e-3)   # logsumexp smoothing for hard constraints

# Edit these three paths to match your local setup:
TROTS_DIR    = Path(globals().get('TROTS_DIR',    r'path/to/TROTS/Protons'))
MATLAB_W_DIR = Path(globals().get('MATLAB_W_DIR', r'path/to/matlab_weights'))
OUTPUT_DIR   = Path(globals().get('OUTPUT_DIR',   r'path/to/output'))

CASE_LABEL = f'TROTS Proton {CASE:02d}'
MAT_FILE   = TROTS_DIR   / f'Protons_{CASE:02d}.mat'
MATLAB_W   = MATLAB_W_DIR / f'Matlab_Proton{CASE:02d}_w.txt'
OUT_FILE   = OUTPUT_DIR  / f'dvh_comparison_{CASE_LABEL.replace(" ", "")}.png'
TXT_FILE   = OUT_FILE.with_suffix('.txt')

_summary = []

# ─── Helpers ─────────────────────────────────────────────────────────────────

def cpp_objective(w, terms):
    """f(w) = sum_i weight_i * sum(max(+/-(D_i*w - bound_i), 0)^2) / nVox
    For MeanDose entries (matRad linear default): f = weight * |mean(D*w) - bound|"""
    total = 0.0
    for D, obj_type, bound, weight in terms:
        dose = D @ w
        n    = len(dose)
        if obj_type == 'MeanDose':
            total += weight * abs(float(dose.mean()) - bound)  # matRad: linear |mean - dRef|
            continue
        if obj_type == 'SquaredOverdosing':
            diff = np.maximum(dose - bound, 0.0)
        elif obj_type == 'SquaredUnderdosing':
            diff = np.maximum(bound - dose, 0.0)
        elif obj_type == 'SquaredDeviation':
            diff = dose - bound
        else:
            continue
        total += weight * float(np.dot(diff, diff)) / n
    return total


def constraint_violations(w, con_terms, eps=CONSTRAINT_EPS):
    """Return list of (con_type, bound, viol_Gy, name) using logsumexp smooth approx.
    Uses the same formula IPOPT optimizes, so a converged solution gives ~0 violation."""
    results = []
    for D, con_type, bound, name in con_terms:
        dose = D @ w
        if con_type == 'MaxDose':
            # smooth max: eps * log(sum(exp((dose - dose.max()) / eps))) + dose.max()
            shift = float(dose.max())
            smooth_max = eps * float(np.log(np.sum(np.exp((dose - shift) / eps)))) + shift
            viol = max(0.0, smooth_max - bound)
        else:  # MinDose: smooth min = -smooth_max(-dose)
            shift = float(dose.min())
            smooth_min = -(eps * float(np.log(np.sum(np.exp((-dose + shift) / eps)))) - shift)
            viol = max(0.0, bound - smooth_min)
        results.append((con_type, bound, viol, name))
    return results

# ─── 1a. Load MATLAB weights (optional) ───────────────────────────────────────

has_matlab = MATLAB_W.exists()
if has_matlab:
    print("Loading MATLAB weights...")
    w_matlab = np.loadtxt(MATLAB_W)
    if w_matlab.size == 0:
        print(f"  MATLAB weight file is empty ({MATLAB_W}) — skipping MATLAB comparison")
        has_matlab = False
        w_matlab   = None
    else:
        print(f"  {len(w_matlab)} bixels  |  range [{w_matlab.min():.3f}, {w_matlab.max():.3f}]")
else:
    print(f"  No MATLAB weights at {MATLAB_W} — skipping MATLAB comparison")
    w_matlab = None

# ─── 1b. Load TROTS reference solution ───────────────────────────────────────

def load_trots_solution(mat_file):
    with h5py.File(mat_file) as f:
        sx = f['solutionX'][()].flatten()
    return sx.astype(np.float64)

print("Loading TROTS reference solution (solutionX)...")
w_trots = load_trots_solution(MAT_FILE)
print(f"  {len(w_trots)} bixels  |  range [{w_trots.min():.3f}, {w_trots.max():.3f}]")

# ─── 2. Load TROTS data ───────────────────────────────────────────────────────

def load_trots(mat_file):
    with h5py.File(mat_file) as f:
        def rstr(ref):
            return ''.join(chr(int(c)) for c in f[ref][()].flatten())

        # Load all D matrices and build name → first-match index (matRad style:
        # readTROTSPatient.m uses find(strcmpi(doseNames, structName), 1) — first match wins).
        d_mats = {}
        name_to_first_id = {}   # lowercase name → 1-based data matrix index
        for idx in range(len(f['data']['matrix']['A'])):
            ref  = f['data']['matrix']['A'][idx, 0]
            mat  = f[ref]
            name = rstr(f['data']['matrix']['Name'][idx, 0]).lower()
            if name not in name_to_first_id:
                name_to_first_id[name] = idx + 1
            if mat.attrs.get('MATLAB_class') == b'double':
                data  = mat['data'][()]
                ir    = mat['ir'][()]
                jc    = mat['jc'][()]
                n_vox = int(mat.attrs['MATLAB_sparse'])
                n_bix = jc.size - 1
                col   = np.zeros(jc[-1], dtype=np.int64)
                for j in range(n_bix):
                    col[jc[j]:jc[j+1]] = j
                D = sps.csr_matrix((data, (ir, col)), shape=(n_vox, n_bix))
            else:
                D = sps.csr_matrix(f[ref][()].T)
            d_mats[idx + 1] = D

        # Replicate matRad's scenario truncation (readTROTSPatient.m lines 320-323):
        # if rows == 9 * svoxels, the matrix contains 9 robustness scenarios — keep only
        # the first svoxels rows (scenario 1 = nominal plan).
        # Read patient.SampledVoxels to get the exact sampled voxel count per structure.
        svoxels_by_name = {}
        n_structs = len(f['patient']['StructureNames'])
        for si in range(n_structs):
            sname = rstr(f['patient']['StructureNames'][si, 0]).lower()
            sv    = f[f['patient']['SampledVoxels'][si, 0]][()]
            svoxels_by_name[sname] = sv.shape[0]  # (N,3) in h5py → N voxels

        print(f"  svoxels: { {k: v for k, v in svoxels_by_name.items() if 'ctv' in k or 'spinal' in k or 'brain' in k} }")
        for idx in list(d_mats.keys()):
            D    = d_mats[idx]
            nm   = rstr(f['data']['matrix']['Name'][idx - 1, 0]).lower()
            svc  = svoxels_by_name.get(nm)
            if svc is not None and D.shape[0] % 9 == 0 and D.shape[0] // 9 == svc:
                print(f"  Truncating [{idx}] '{nm}': {D.shape[0]} → {svc} rows")
                d_mats[idx] = D[:svc, :]

        # Build reverse map: 1-based index → lowercase matrix name
        idx_to_matname = {}
        for idx in range(len(f['data']['matrix']['A'])):
            idx_to_matname[idx + 1] = rstr(f['data']['matrix']['Name'][idx, 0]).lower()

        # matRad only creates objectives for entries whose name matches a patient structure.
        # Non-structure entries (regularization, Total MU, etc.) are silently skipped.
        structure_names_lc = {rstr(f['patient']['StructureNames'][si, 0]).lower()
                              for si in range(n_structs)}

        entries = []
        for i in range(len(f['problem']['dataID'])):
            name     = rstr(f['problem']['Name'][i, 0])
            if name.lower() not in structure_names_lc:
                continue   # non-structure entry — matRad skips these
            orig_id  = int(f[f['problem']['dataID'][i, 0]][()].flatten()[0])
            orig_mat_name = idx_to_matname.get(orig_id, '')
            # matRad detects MeanDose objectives by checking if data.matrix(dataID).Name
            # equals structure_name + " (mean)".  Use the regular (non-mean) D matrix
            # for dose computation in both cases — same as matRad's dij.physicalDose.
            is_mean = (orig_mat_name == name.lower() + ' (mean)')
            # Use the same D matrix matRad uses: first name-matched, not dataID-matched.
            data_id  = name_to_first_id.get(name.lower(), orig_id)
            if data_id != orig_id:
                D_matched = d_mats[data_id]
                D_orig    = d_mats[orig_id]
                if (D_matched.shape[0] % 9 == 0 and
                        D_matched.shape[0] // 9 == D_orig.shape[0]):
                    d_mats[data_id] = D_matched[:D_orig.shape[0], :]
            entries.append({
                'name':          name,
                'data_id':       data_id,
                'is_constraint': bool(f[f['problem']['IsConstraint'][i, 0]][()].flatten()[0]),
                'minimize':      bool(f[f['problem']['Minimise'][i, 0]][()].flatten()[0]),
                'weight':        float(f[f['problem']['Weight'][i, 0]][()].flatten()[0]),
                'bound':         float(f[f['problem']['Objective'][i, 0]][()].flatten()[0]),
                'is_mean':       is_mean,
            })

        return entries, d_mats


print(f"\nLoading {CASE_LABEL}...")
entries, d_mats = load_trots(MAT_FILE)
n_bixels = d_mats[1].shape[1]
print(f"  {len(entries)} problem entries, {n_bixels} bixels")

if has_matlab:
    assert len(w_matlab) == n_bixels, \
        f"MATLAB weight vector length {len(w_matlab)} != {n_bixels} bixels — wrong CASE?"

# ─── 3. Build optimization problem ───────────────────────────────────────────

# TROTS semantics:
#   is_constraint=False, minimize=True  → OAR soft obj → SquaredOverdosing(bound=d_max)
#   is_constraint=False, minimize=False → Target soft obj → SquaredUnderdosing(bound=d_min)
#   is_constraint=True,  minimize=True  → hard max dose constraint (keep dose ≤ bound)
#   is_constraint=True,  minimize=False → hard min dose constraint (require dose ≥ bound)

struct_terms = []    # (D, obj_type, bound, weight)   obj_type ∈ {SquaredOverdosing, SquaredUnderdosing, MeanDose}
constraint_terms = []  # (D, con_type, bound, name)
mean_entries = []    # names of mean-dose soft objectives (for logging)

for e in entries:
    D  = d_mats[e['data_id']]

    if e['is_constraint']:
        con_type = 'MaxDose' if e['minimize'] else 'MinDose'
        constraint_terms.append((D, con_type, e['bound'], e['name']))
    else:
        if e['is_mean']:
            # matRad uses matRad_MeanDose: penalize mean(D@w) > bound
            obj_type = 'MeanDose'
            mean_entries.append(e['name'])
        else:
            obj_type = 'SquaredOverdosing' if e['minimize'] else 'SquaredUnderdosing'
        struct_terms.append((D, obj_type, e['bound'], e['weight']))

print(f"  {len(struct_terms)} soft objectives, {len(constraint_terms)} hard constraints")
if mean_entries:
    print(f"  MeanDose objectives: {mean_entries}")

# Build name → {col: [bounds]} for each objective/constraint type
_COL_ORDER = ['OD obj', 'UD obj', 'Mean obj', 'MaxDose', 'MinDose']
_name_bounds = {}
for e in entries:
    if e['is_constraint']:
        col = 'MaxDose' if e['minimize'] else 'MinDose'
    elif e['is_mean']:
        col = 'Mean obj'
    else:
        col = 'OD obj' if e['minimize'] else 'UD obj'
    _name_bounds.setdefault(e['name'], {}).setdefault(col, []).append(e['bound'])
_used_cols = [c for c in _COL_ORDER if any(c in v for v in _name_bounds.values())]

x0 = np.zeros(n_bixels)  # values ignored — C++ auto-scales from MinDose bounds

f_trots_val = cpp_objective(w_trots, struct_terms)
_trots_line  = f" Reference |  f={f_trots_val:.6f}"
print("TROTS reference:     " + _trots_line)
_summary.append("TROTS reference:     " + _trots_line)

if has_matlab:
    f_mat_val   = cpp_objective(w_matlab, struct_terms)
    _matlab_line = f" Converged |  f={f_mat_val:.6f}"
    print("MATLAB IPOPT (matRad):" + _matlab_line)
    _summary.append("MATLAB IPOPT (matRad):" + _matlab_line)

    print("  Objective breakdown at MATLAB weights (top contributors):")
    _term_vals = []
    for D, obj_type, bound, weight in struct_terms:
        dose = D @ w_matlab
        if obj_type == 'MeanDose':
            v = weight * abs(float(dose.mean()) - bound)
        elif obj_type == 'SquaredOverdosing':
            diff = np.maximum(dose - bound, 0.0)
            v = weight * float(np.dot(diff, diff)) / len(dose)
        else:
            diff = np.maximum(bound - dose, 0.0)
            v = weight * float(np.dot(diff, diff)) / len(dose)
        _term_vals.append((v, obj_type, bound, weight, D.shape[0]))
    for v, ot, b, w_v, n in sorted(_term_vals, reverse=True)[:10]:
        print(f"    {ot:22s} bound={b:6.1f} Gy  n={n:5d}  weight={w_v:.4g}  f_term={v:.1f}")

# ─── 4. Run C++ IPOPT ────────────────────────────────────────────────────────

print(f"\nRunning C++ IPOPT optimizer (max_iter={MAX_ITER})...")
try:
    import PythonQt
    our_opt = PythonQt.qSlicerExternalBeamPlanningModuleWidgets.qSlicerIpoptOptimizer()
    our_opt.setMaxIterations(MAX_ITER)
    our_opt.setOption("limited_memory_max_history", LBFGS_HISTORY)
    our_opt.setOption("acceptable_tol",             ACCEPTABLE_TOL)
    our_opt.setOption("acceptable_obj_change_tol",  ACCEPTABLE_OBJ_CHANGE)
    our_opt.setOption("linear_solver",              LINEAR_SOLVER)
    our_opt.clearStructureTerms()

    for D, obj_type, bound, weight in struct_terms:
        D_coo = D.tocoo()
        our_opt.addStructureTerm(
            D_coo.data.tolist(),
            D_coo.row.astype(int).tolist(),
            D_coo.col.astype(int).tolist(),
            int(D.shape[0]), int(D.shape[1]),
            obj_type, float(bound), float(weight)
        )

    for D, con_type, bound, _cname in constraint_terms:
        D_coo = D.tocoo()
        our_opt.addStructureConstraint(
            D_coo.data.tolist(),
            D_coo.row.astype(int).tolist(),
            D_coo.col.astype(int).tolist(),
            int(D.shape[0]), int(D.shape[1]),
            con_type, float(bound), CONSTRAINT_EPS
        )

    t0         = time.time()
    success    = our_opt.solve(x0.tolist())
    t_ours     = time.time() - t0
    w_ours     = np.asarray(our_opt.getSolution())
    status_str = 'Converged' if success else 'MAX_ITER reached / not converged'
    has_ours   = True
except Exception as exc:
    print(f"  qSlicerIpoptOptimizer unavailable ({exc})")
    has_ours   = False
    t_ours     = 0.0
    w_ours     = x0.copy()

f_ours_val = cpp_objective(w_ours, struct_terms)
if has_ours:
    _line = f" {status_str}  |  f={f_ours_val:.6f}  |  t={t_ours:.1f}s"
    print("C++ IPOPT (SlicerRT):" + _line)
    _summary.append("C++ IPOPT (SlicerRT):" + _line)

# ─── 5. Collect structures for DVH ───────────────────────────────────────────

DVH_STRUCTURES = {}
for e in entries:
    if e['name'] not in DVH_STRUCTURES and d_mats[e['data_id']].shape[0] >= 100:
        DVH_STRUCTURES[e['name']] = e['data_id']

# ─── 6. Compute and plot DVH ──────────────────────────────────────────────────

n_structs = len(DVH_STRUCTURES)
n_cols    = 4
n_rows    = (n_structs + 1 + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4.5, n_rows * 4.5))
axes = axes.flatten()

label_ours   = 'C++ IPOPT (SlicerRT)' if has_ours else 'C++ IPOPT (unavailable)'
label_matlab = 'MATLAB IPOPT (matRad)'

title_vs = f' vs {label_matlab}' if has_matlab else ''
fig.suptitle(f'DVH Comparison: {label_ours}{title_vs}\n({CASE_LABEL})', fontsize=13)

for ax, (name, did) in zip(axes, DVH_STRUCTURES.items()):
    D       = d_mats[did]
    d_cpp   = D @ w_ours
    d_trots = D @ w_trots

    vmax = max(d_cpp.max(), d_trots.max())
    if has_matlab:
        d_mat = D @ w_matlab
        vmax  = max(vmax, d_mat.max())
    vmax = vmax * 1.05
    bins = np.linspace(0, vmax, 300)

    ax.plot(bins, np.array([np.mean(d_cpp   >= b) for b in bins]) * 100, 'r-',  lw=2, label=label_ours)
    if has_matlab:
        ax.plot(bins, np.array([np.mean(d_mat   >= b) for b in bins]) * 100, 'b--', lw=2, label=label_matlab)
    ax.plot(bins, np.array([np.mean(d_trots >= b) for b in bins]) * 100, 'g:',  lw=2, label='TROTS reference')
    ax.set_title(name, fontsize=9)
    ax.set_xlabel('Dose (Gy)', fontsize=8)
    ax.set_ylabel('Volume (%)', fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 105)

for ax in axes[n_structs:]:
    ax.set_visible(False)

ax_sum = axes[n_structs]
ax_sum.set_visible(True)
ax_sum.axis('off')
sum_rows = [
    [label_ours,        f'{f_ours_val:.4f}'],
    ['TROTS reference', f'{f_trots_val:.4f}'],
]
if has_matlab:
    sum_rows.insert(1, [label_matlab, f'{f_mat_val:.4f}'])
tbl = ax_sum.table(cellText=sum_rows, colLabels=['Solver', 'Objective f(w)'],
                   loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 2)
ax_sum.set_title('Summary', fontsize=10)

plt.tight_layout()
plt.savefig(OUT_FILE, dpi=150, bbox_inches='tight')
print(f"\nDVH comparison saved to: {OUT_FILE}")

# ─── 7. Text summary ─────────────────────────────────────────────────────────

_summary.append(f"\nObjective  f(w) = sum_i penalty_i * sum(max(+/-diff,0)^2)/n  [C++ formula]:")
_summary.append(f"  {label_ours:30s}  f = {f_ours_val:.6f}")
if has_matlab:
    _summary.append(f"  {label_matlab:30s}  f = {f_mat_val:.6f}")
_summary.append(f"  {'TROTS reference':30s}  f = {f_trots_val:.6f}")

if constraint_terms:
    viol_ours  = constraint_violations(w_ours,  constraint_terms)
    viol_trots = constraint_violations(w_trots, constraint_terms)
    viol_pairs = [(label_ours, viol_ours)]
    if has_matlab:
        viol_pairs.insert(1, (label_matlab, constraint_violations(w_matlab, constraint_terms)))
    viol_pairs.append(('TROTS reference', viol_trots))
    _summary.append(f"\nHard constraint violations (logsumexp, >0.01 Gy threshold):")
    for label, viol_list in viol_pairs:
        n_v = sum(1 for _, _, v, _ in viol_list if v > 0.01)
        mv  = max(v for _, _, v, _ in viol_list)
        _summary.append(f"  {label:30s}  {n_v:2d} violated  |  max = {mv:.3f} Gy")
        for con_type, bound, viol, name in sorted(viol_list, key=lambda x: -x[2]):
            if viol > 0.01:
                _summary.append(f"    {label[:6]:6s}  {con_type:8s} {name:35s} bound={bound:.1f} Gy  viol={viol:.3f} Gy")

_col_w = 10
_matlab_col = f"  {'MATLAB':>10s}" if has_matlab else ''
_header = (f"{'Structure':35s}  {'C++ IPOPT':>10s}{_matlab_col}  {'TROTS ref':>10s}"
           + ''.join(f"  {c:>{_col_w}s}" for c in _used_cols))
_summary.append(f"\nMean dose summary (Gy):")
_summary.append(_header)
_summary.append("-" * len(_header))
for name, did in DVH_STRUCTURES.items():
    D        = d_mats[did]
    md_cpp   = float((D @ w_ours).mean())
    md_trots = float((D @ w_trots).mean())
    bounds = _name_bounds.get(name, {})
    cols   = ''.join(
        f"  {', '.join(f'{b:.1f}' for b in bounds[c]):>{_col_w}s}" if c in bounds else f"  {'---':>{_col_w}s}"
        for c in _used_cols)
    matlab_col = f'  {float((D @ w_matlab).mean()):10.3f}' if has_matlab else ''
    _summary.append(
        f"{name:35s}  {md_cpp:10.3f}{matlab_col}  {md_trots:10.3f}{cols}")

print('\n'.join(_summary))
with open(TXT_FILE, 'w') as fh:
    fh.write('\n'.join(_summary) + '\n')
print(f"Summary saved to: {TXT_FILE}")
