PY ?= python

.PHONY: install test figures

install:
	$(PY) -m pip install -e ".[dev]"

test:
	$(PY) -m pytest -q

figures:
	$(PY) scripts/phase0_benchmark.py
