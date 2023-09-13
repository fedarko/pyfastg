.PHONY: test stylecheck style

test:
	python3 -B -m pytest pyfastg/tests/ --cov-report xml --cov-report term

stylecheck:
	flake8 pyfastg setup.py
	black --check -l 79 pyfastg setup.py

style:
	black -l 79 pyfastg setup.py
