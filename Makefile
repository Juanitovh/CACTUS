.PHONY: all build-cpp install install-dev clean

all: build-cpp install

build-cpp:
	$(MAKE) -C source_fibre_optimization

install: build-cpp
	pip install -e .

install-dev: build-cpp
	pip install -e ".[dev]"

clean:
	$(MAKE) -C source_fibre_optimization clean
