.PHONY: all test dist

all: test

test:
	pytest

dist:
	python -m build
