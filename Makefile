PYTHON ?= python3

.PHONY: install test test-fast validate figures sweeps notebooks report all clean

install:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	$(PYTHON) -m pytest

test-fast:
	$(PYTHON) -m pytest -m "not slow"

validate:
	$(PYTHON) scripts/run_validation_suite.py

figures:
	$(PYTHON) scripts/generate_all_figures.py

sweeps:
	$(PYTHON) scripts/run_parameter_sweeps.py

notebooks:
	$(PYTHON) scripts/build_notebooks.py

report:
	$(PYTHON) scripts/reproduce_report.py

all: validate sweeps figures report

clean:
	rm -rf .pytest_cache build dist src/*.egg-info
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
