.PHONY: stylecheck style

stylecheck:
	flake8 pyfastg
	black --check -l 79 pyfastg

style:
	black -l 79 pyfastg
