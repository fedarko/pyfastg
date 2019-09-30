.PHONY: test stylecheck style

test:
	python3 -B -m pytest pyfastg/tests/ --cov

stylecheck:
	flake8 pyfastg
	black --check -l 79 pyfastg

style:
	black -l 79 pyfastg
