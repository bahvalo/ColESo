CXX?=c++
CCX?=cc
FCX?=gfortran
COMPILER_OPTIONS = -O3
INCLUDES = -I../src
OUTDIR = ../bin

COLESO_LIB = $(OUTDIR)/coleso.a
COLESO_LIBQD = $(OUTDIR)/colesoqd.a

testc = $(OUTDIR)/testc
testcpp = $(OUTDIR)/testcpp
testcppqd = $(OUTDIR)/testcppqd
testf = $(OUTDIR)/testf

testc: $(testc);
testcpp: $(testcpp);
testcppqd: $(testcppqd);
testf: $(testf);

default: testcpprun ;

$(testc): test.c $(COLESO_LIB)
	@echo  compiling $@ ...
	@$(CCX) $(COMPILER_OPTIONS) $(INCLUDES) test.c $(COLESO_LIB) -lstdc++ -lm -o $@

$(testcpp): test.cpp $(COLESO_LIB)
	@echo  compiling $@ ...
	@$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) test.cpp $(COLESO_LIB) -o $@

$(testcppqd): test.cpp $(COLESO_LIBQD)
	@echo  compiling $@ ...
	@$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) -I.. -DEXTRAPRECISION_COLESO test.cpp $(COLESO_LIBQD) $(OUTDIR)/libqd.a -o $@

$(testf): test.f $(COLESO_LIB)
	@echo  compiling $@ ...
	@$(FCX) test.f $(INCLUDES) $(COLESO_LIB) -J $(OUTDIR) -lstdc++ -o $(OUTDIR)/$@

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
	rm -rf $(testc) $(testcpp) $(testcppqd) $(testf) ./DATA1D ./DATA2D

