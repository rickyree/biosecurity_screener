# Session Setup

**CRITICAL: Work in the user's current directory.** The working directory where Claude Code was invoked is the project root. Never create files in `~/` or other arbitrary locations.

## 1. Understand Project Context First

Before creating any files or directories, examine what already exists:

```bash
pwd                             # Confirm working directory
ls -la                          # See existing structure
```

Look for: existing project organization, input files (PDBs, FASTAs), config files, previous Amina outputs. Adapt to the existing structure rather than imposing a new one.

## 2. Ensure Amina is Up to Date

```bash
pip install --upgrade amina-cli           # or: uv pip install --upgrade amina-cli
```

This is a no-op if already current.

## 3. Dependency Management

**All packages must be installed in the same environment as amina-cli (Python 3.11+).**

```bash
# Verify environment
python --version                          # Must be 3.11+
which python                              # Confirm correct interpreter
```

**Install dependencies using the same method as amina-cli:**

```bash
# If amina was installed with pip:
pip install biopython matplotlib

# If amina was installed with uv:
uv pip install biopython matplotlib
```

If Python is <3.11 or amina is not found, the user needs to activate the correct environment or install amina-cli first:

```bash
# Example setup (if needed)
uv venv --python 3.11                     # or: python3.11 -m venv .venv
source .venv/bin/activate
pip install amina-cli                     # or: uv pip install amina-cli
```

## 4. Verify Amina Authentication (silently)

Only prompt user if something fails:

```bash
amina auth status               # If unauthenticated: amina auth set-key "{key}"
```

## 5. Output Session Summary

Detect the active environment:

```bash
echo $VIRTUAL_ENV                         # Set by venv/uv venv when activated
which python                              # Shows interpreter path
```

- **venv/uv venv**: `$VIRTUAL_ENV` is set (e.g., `/path/to/project/.venv`)
- **uv run** (ephemeral): `$VIRTUAL_ENV` empty but python points to `.venv` or uv cache
- **System Python**: `$VIRTUAL_ENV` empty, python points to `/usr/bin/python` or similar

After completing setup checks, output a brief summary:

```
**Amina Session Ready**
- **Directory**: {cwd}
- **Environment**: {$VIRTUAL_ENV path, "uv-managed", or "system Python"}
- **Python**: {version} | **Amina**: {version} ✓
- **Auth**: {authenticated/needs setup}
- **Project files**: {e.g., "3 PDB files, 1 FASTA" or "empty directory"}
```

Keep it concise. Only mention issues if action is needed.
