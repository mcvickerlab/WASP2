# Regression tests

Run `pytest tests/regression/ -v` for the mapping regressions. The deterministic
SNV and indel fixtures live alongside the tests and require the native extension.
The installed CLI acceptance test also requires samtools, bcftools, bedtools and
BWA. It checks exact counts and swapped sequences, remapping/filtering, analysis
identities, and single-cell matrices through `scripts/smoke_test.sh`.

Historical benchmark checks in `test_pipeline_regression.py` require optional
local `baselines/` measurements and skip explicitly when those are absent. Their
MD5 and timing checks are historical evidence, not the release acceptance gate.

Bulk and single-cell statistical reference comparisons are under `tests/analysis`
and `tests/counting`. They execute the pinned Python archive and distinguish its
narrow-integer limitations from explicitly widened reference inputs. Keep exact
count comparisons separate from statistical tolerances.
