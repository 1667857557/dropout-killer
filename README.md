# DropoutKiller

Selective dropout recovery for single-cell RNA-seq.

## Design

The package follows a zero-preserving strategy:

1. Build biological local neighborhoods using SuperCell-style membership.
2. Estimate dropout confidence using local low-rank reconstruction.
3. Modify only high-confidence technical dropout events.
4. Recover expression using local cell borrowing and optional gene-network priors.

## Validation

Tests cover:

- matrix dimension preservation
- selective recovery output integrity
- synthetic dropout workflow

## Status

Development version. Large-scale sparse optimization and full weighted-neighbor implementation are ongoing.
