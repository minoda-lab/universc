SHELL := /bin/bash -O extglob

VERSION=$(shell git describe --tags --always --dirty)
export ROOT_DIR=$(shell pwd)

# Targets for development builds.
#
all: reference

clean: reference-clean

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
	echo "removing pre-generated cellranger reference files"
	rm -Rf  test/cellranger_reference/cellranger-tiny-ref/3.0.0
	rm -Rf  test/cellranger_reference/cellranger-tiny-ref/1.2.0

manual:
	# add manual directory to PATH if not already found
	## check config for Linux
	if [[ -f /etc/manpath.config ]]; then CONFIG=$(echo /etc/manpath.config); fi
	## check config for Mac
	if [[ -f /etc/manpaths ]]; then CONFIG:=$(echo /etc/manpaths); fi
	echo $(CONFIG)
	#if [[ ! -z $CONFIG ]];  then MANDIR=`tail -n 1 ${CONFIG}`; else if [[ ! -z $MANPATH ]]; then SHELL_RC=`echo ~/.${0}rc`; echo "export MANPATH=/usr/local/man" >> $SHELL_RC; MANDIR=`echo ${MANPATH} | cut -d: -f1`; fi; fi
	#mkdir -p ${MANDIR}/man1
	#cp man/convert.sh man/convert.sh.1
	#mv man/convert.sh.1 ${MANDIR}/man1

CONFIG:="/etc/manpaths"

.PHONY: print_vars

print_vars:
	echo $(CONFIG)	
