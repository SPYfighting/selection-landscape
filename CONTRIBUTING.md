# Contributing

Open an issue with a minimal input, the command used, the expected result,
the observed result, and the software versions. Use synthetic examples when
reporting unpublished data or sequences.

Install with `python -m pip install -e '.[test]'` and run `python -m pytest`.
Changes to filtering, score definitions, calibration, candidate eligibility,
or ordering require a test against an independently specified expected result.
Numerical changes should be described separately from file-format or interface
changes. Keep score construction independent of model fitting and prediction.

Do not commit raw reads, study labels, fitted study models, credentials, or
machine-specific paths. Each contribution must retain applicable third-party
copyright and license notices.
