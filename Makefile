SHELL := /bin/bash -O extglob

VERSION=$(shell git describe --tags --always --dirty)
export ROOT_DIR=$(shell pwd)

# Targets for development builds.
#
all: reference

clean: reference-clean manual-clean

# Targets for python builds

# Find python libs, includes, and default flags.
PYTHON_ROOT:=$(shell python-config --exec-prefix)
PYTHON_LIB_DIR:=$(dir $(lastword $(wildcard $(PYTHON_ROOT)/lib*/libpython*.so*)))
PYTHON_SITE_PKG:=$(lastword $(wildcard $(PYTHON_LIB_DIR)/python*/site-packages))
PYTHON_CFLAGS:=$(shell python-config --cflags)
PYTHON_LDFLAGS:=$(shell python-config --ldflags)

reference:
	make -C  test/cellranger_reference/cellranger-tiny-ref

reference-clean:
	make -C  test/cellranger_reference/cellranger-tiny-ref clean
	@echo "removing pre-generated cellranger reference files"
	rm -Rf  test/cellranger_reference/cellranger-tiny-ref/3.0.0
	rm -Rf  test/cellranger_reference/cellranger-tiny-ref/1.2.0

# Copy to manuals (requires root priviledges
manual:
	@bash man/INSTALL

manual-clean:
	@bash man/REMOVE

install:
	@bash inst/INSTALL prefix=$(prefix)

remove:
	@bash inst/REMOVE

uninstall: remove clean

reinstall: uninstall install
