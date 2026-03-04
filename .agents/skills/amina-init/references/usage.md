# Usage Guide

## Agent Behavior

**Hands-free execution**: Run all commands automatically without asking for permission. Never ask "Want me to...?" or "Should I...?" — just do it. Use defaults and infer parameters from context. The only exceptions are:
- Truly ambiguous requests where you cannot infer the user's intent
- Destructive operations (deleting files, overwriting important data)
- Missing required inputs that cannot be inferred (e.g., target protein for binder design)

**Complete workflows end-to-end**: When a job is submitted, automatically:
1. Monitor its status (poll periodically or use `amina jobs wait`)
2. Download results when complete
3. Parse and present the outputs
4. Suggest next steps

Do NOT stop mid-workflow to ask if the user wants to continue. Assume yes unless the task is complete or an error occurs.

**Adapt to expertise**: Jargon (scRMSD, pLDDT, contigs) → be terse. Plain English → explain steps.

**Clarity in communication**: Output information in a way that LLMs can easily parse and understand. For user facing information, use markdown formatting and emojis sparingly.

**Parse requests for**: Target (PDB ID/sequence/file), Goal, Constraints (hotspots, chains, lengths), Scale (number of designs).

**Classify workflow**:
- "binder", "bind to" → De Novo Binder Design
- "design protein" → De Novo Monomer Design
- "motif", "scaffold", "contigs" → Motif Scaffolding
- "enzyme", "catalytic", "mutant" → Enzyme Engineering
- "dock", "ligand", "SMILES" → Small Molecule Docking
- "validate", "assess" → Structure Validation
- "interface", "PPI" → PPI Analysis
- "embedding", "cluster" → Sequence Space Exploration

**Research when useful**: For unfamiliar proteins, use `WebSearch` for PDB IDs, binding interfaces, known hotspots. Skip for self-contained requests.

## CLI Reference

```bash
# Discovery (always check before running)
amina tools list                # List all tools
amina run <tool> --help         # Show parameters (never guess)

# Execution
amina run <tool> [params] -o <dir>           # Sync
amina run <tool> --background --job-name <n> # Async (>30s jobs)

# Job management
amina jobs list --status running
amina jobs status <job-id>
amina jobs wait <job-id> && amina jobs download <job-id> -o ./
amina jobs logs <job-id>        # Debug failures
```

**Use `--background` for**: Structure prediction (ESMFold, Boltz), RFDiffusion, MD simulations, any job >30s.

## Execution Patterns

### Project Organization

**Always use relative paths from the current working directory.** Never use absolute paths like `~/` or `/Users/...`.

```bash
# All outputs go in current directory or subdirectories
./                              # Current working directory is project root
./downloads/                    # Downloaded files (e.g., PDBs, FASTAs)
./results/<tool>_<jobname>/     # Single tool runs
./01_structure/ 02_design/      # Multi-stage workflows
./round_1/ round_2/             # Iterative refinement
```

If the directory already has structure, follow it. If empty, create minimal organization as needed.

### Parallelization

Use `--background` for ALL parallel jobs, both long and short:

1. **Submit**: Launch all jobs with `--background` in a single message (multiple Bash tool calls — each returns a job ID instantly)
2. **Wait**: `amina jobs wait {id1} {id2} {id3} ... --poll-interval 10` (blocks until all complete)
3. **Download**: `amina jobs download {id} -o {dir}/` for each job

### Tool Chaining

| From | Output | To | Input |
|------|--------|-----|-------|
| RFdiffusion | `*.pdb` | ProteinMPNN | `--pdb` |
| ProteinMPNN | `*.fasta` | ESMFold | `--sequence` (parse FASTA) |
| ESMFold | `*.pdb` | Simple RMSD | `--mobile` |
| P2Rank | `*_predictions.csv` | Vina | `--center-x/y/z` (parse CSV) |
| PDB Cleaner | `*_cleaned.pdb` | Any PDB tool | `--pdb` |

After `amina run`, use `Glob` to find outputs. Never guess metrics—wait for actual results.

## Results & Iteration

Present results naturally for the task. Read output files and interpret for user. Suggest 1-2 follow-ups conversationally.

**Key principle**: Never rerun expensive stages. Check what exists before re-executing.

## Claude Skills

You have access to many claude skills that you will need for protein engineering tasks. These skills are installed in the `.claude/skills/` directory.

There are 3 types of skills:
- **Amina Tool Skills**: Skills for using the Amina tools e.g. amina-rfdiffusion, amina-proteinmpnn, amina-boltz-2, amina-diffdock (prefixed with "amina-")
- **Workflow Skills**: Skills for running end-to-end workflows e.g., workflow-binder-design, workflow-enzyme-engineering (prefixed with "workflow-")
- **External Tool Skills**: Skills for using external tools e.g., pdb-database, uniprot-database, pymol, biopython, etc. (no prefix)

These skills should always be used before running any tool or workflow to ensure you are using the tools with scientific best practices.

## Error Handling

```bash
amina auth status                                # Auth issues
amina jobs logs <job-id-1> <job-id-2> ...        # Job failures
amina tools <tool> --help                        # Format/parameter issues
```
