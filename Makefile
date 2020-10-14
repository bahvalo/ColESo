CXX?=c++
CCX?=cc
FCX?=gfortran
COMPILER_OPTIONS = -O3
INCLUDES = -I./src
OUTDIR = ./bin
QDSUFF ?=
LIBS ?= $(OUTDIR)/coleso.a

default: examplecpp ;

qd-package:
	$(MAKE) -C ./qd

lib:
	$(MAKE) -C ./src

libqd:
	$(MAKE) -C ./src CFLAGS=-DEXTRAPRECISION_COLESO QDSUFF=qd

examplecpp: lib
	$(MAKE) -C ./example examplecpp

examplecppqd: libqd qd-package
	$(MAKE) -C ./example examplecppqd

examplec: lib
	$(MAKE) -C ./example examplec

examplef: lib
	$(MAKE) -C ./example examplef

testcpp: lib
	$(MAKE) -C ./test testcpp

testc: lib
	$(MAKE) -C ./test testc

testf: lib
	$(MAKE) -C ./test testf

testcpprun: lib
	$(MAKE) -C ./test testcpprun

testcrun: lib
	$(MAKE) -C ./test testcrun

testfrun: lib
	$(MAKE) -C ./test testfrun

check: testcpp
	$(OUTDIR)/testcpp test

checkqd: libqd qd-package
	$(MAKE) -C ./test testcppqd
	$(OUTDIR)/testcppqd test waveinchannel-succ
	$(OUTDIR)/testcppqd test waveinchannel-succ-dd
	$(OUTDIR)/testcppqd test waveinchannel-succ-qd

clean:
	rm -rf bin
	$(MAKE) -C ./test clean

