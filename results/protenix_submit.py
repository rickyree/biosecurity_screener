#!/usr/bin/env python3
"""
Submit one Boltz2 job per row in candidate_sequences.tsv,
collect job IDs, wait, download, then run dMaSIF workflow on results.
"""

import csv, os, re, subprocess, sys, glob, shutil

TSV    = 'results/dmasif/candidate_sequences.tsv'
OUTDIR = 'results/boltz2'
os.makedirs(OUTDIR, exist_ok=True)

# ── 1. Read sequences ────────────────────────────────────────────────────────
rows = []
with open(TSV) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        name  = row['candidate']
        seqs  = [row[k] for k in reader.fieldnames[1:] if row[k].strip()]
        if seqs:
            rows.append((name, seqs))

print(f"Found {len(rows)} candidates to submit\n")

# ── 2. Submit all jobs in background ─────────────────────────────────────────
job_ids   = {}   # candidate -> job_id
job_names = {}   # candidate -> amina job_name (sanitised)

for name, seqs in rows:
    # Amina job names: alphanumeric + hyphens only
    jname = re.sub(r'[^A-Za-z0-9\-]', '-', name)[:48]

    cmd = ['amina', 'run', 'boltz2',
           '--job-name', jname,
           '--output', OUTDIR,
           '--background']
    for seq in seqs:
        cmd += ['-s', seq]

    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    # Extract job ID from output (UUID pattern)
    m = re.search(r'Job submitted:\s*([0-9a-f\-]{36})', out)
    if m:
        jid = m.group(1)
        job_ids[name]   = jid
        job_names[name] = jname
        print(f"  submitted  {name:<32}  job={jid[:8]}")
    else:
        print(f"  FAILED     {name:<32}\n  stdout: {result.stdout[:200]}\n  stderr: {result.stderr[:200]}")

print(f"\nSubmitted {len(job_ids)}/{len(rows)} jobs")

# Save job ID map for recovery
with open(f'{OUTDIR}/job_ids.tsv', 'w') as f:
    f.write('candidate\tjob_id\tjob_name\n')
    for name, jid in job_ids.items():
        f.write(f"{name}\t{jid}\t{job_names[name]}\n")
print(f"Job IDs saved → {OUTDIR}/job_ids.tsv")

# ── 3. Wait for all jobs ──────────────────────────────────────────────────────
if not job_ids:
    sys.exit("No jobs submitted.")

all_ids = list(job_ids.values())
print(f"\nWaiting for {len(all_ids)} jobs …")
wait_cmd = ['amina', 'jobs', 'wait', '--poll-interval', '20'] + all_ids
subprocess.run(wait_cmd)

# ── 4. Download results ───────────────────────────────────────────────────────
print("\nDownloading results …")
for name, jid in job_ids.items():
    dl_dir = os.path.join(OUTDIR, name)
    os.makedirs(dl_dir, exist_ok=True)
    r = subprocess.run(['amina', 'jobs', 'download', jid, '-o', dl_dir],
                       capture_output=True, text=True)
    if r.returncode == 0:
        print(f"  ✓  {name}")
    else:
        print(f"  ✗  {name}: {r.stderr[:100]}")

# ── 5. Collect best PDB per candidate → results/boltz2/pdbs/ ─────────────────
PDB_DIR = os.path.join(OUTDIR, 'pdbs')
os.makedirs(PDB_DIR, exist_ok=True)

collected = []
for name, _ in rows:
    dl_dir = os.path.join(OUTDIR, name)
    # Boltz2 names structures like: boltz2_{jobname}_structure.pdb
    pdbs = sorted(glob.glob(os.path.join(dl_dir, '**', 'boltz2_*_structure.pdb'), recursive=True))
    if not pdbs:
        pdbs = sorted(glob.glob(os.path.join(dl_dir, '**', '*.pdb'), recursive=True))
    if pdbs:
        best = pdbs[0]
        dest = os.path.join(PDB_DIR, f"{name}_boltz2.pdb")
        shutil.copy(best, dest)
        collected.append(dest)
        print(f"  ✓  {name} → {os.path.basename(dest)}")
    else:
        print(f"  ✗  {name}: no PDB found in {dl_dir}")

print(f"\nCollected {len(collected)} predicted structures → {PDB_DIR}/")

# ── 6. Run dMaSIF workflow ────────────────────────────────────────────────────
print("\nRunning dMaSIF workflow on predicted structures …\n")
dmasif_cmd = ['python', 'results/dmasif_workflow.py',
              '--use-boltz2',
              '--candidates'] + collected
subprocess.run(dmasif_cmd)
