CXX?=c++
CCX?=cc
FCX?=gfortran
COMPILER_OPTIONS = -O3
INCLUDES = -I../src
OUTDIR = ../bin
LIBS ?= $(OUTDIR)/coleso.a

default: testcpprun ;

testc: $(OUTDIR)/testc
$(OUTDIR)/testc:	
	$(CCX) $(COMPILER_OPTIONS) $(INCLUDES) test.c $(LIBS) -lstdc++ -lm -o $@

testcpp: $(OUTDIR)/testcpp
$(OUTDIR)/testcpp:
	$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) test.cpp $(LIBS) -o $@

testcppqd: $(OUTDIR)/testcppqd
$(OUTDIR)/testcppqd:
	$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) -I.. -DEXTRAPRECISION_COLESO test.cpp $(OUTDIR)/colesoqd.a $(OUTDIR)/libqd.a -o $@

testf: $(OUTDIR)/testf
$(OUTDIR)/testf:
	$(FCX) test.f $(INCLUDES) $(LIBS) -J $(OUTDIR) -lstdc++ -o $(OUTDIR)/$@

testcrun: testc
	@mkdir -p DATA1D DATA2D
	$(OUTDIR)/testc
	sh ./plotall_1d.sh
	sh ./plotall_2d.sh

testcpprun: testcpp
	@mkdir -p DATA1D DATA2D
	$(OUTDIR)/testcpp
	sh ./plotall_1d.sh
	sh ./plotall_2d.sh

testfrun: testf
	@mkdir -p DATA1D DATA2D
	$(OUTDIR)/testf
	sh ./plotall_1d.sh
	sh ./plotall_2d.sh

clean:
	rm -rf $(OUTDIR)/testc $(OUTDIR)/testcpp $(OUTDIR)/testf ./DATA1D ./DATA2D

