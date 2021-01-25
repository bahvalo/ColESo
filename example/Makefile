CXX?=c++
CCX?=cc
FCX?=gfortran
COMPILER_OPTIONS = -O3
INCLUDES = -I../src
OUTDIR = ../bin

COLESO_LIB = $(OUTDIR)/coleso.a
COLESO_LIBQD = $(OUTDIR)/colesoqd.a

examplec = $(OUTDIR)/examplec
examplecpp = $(OUTDIR)/examplecpp
examplecppqd = $(OUTDIR)/examplecppqd
examplef = $(OUTDIR)/examplef

examplec: $(examplec);
examplecpp: $(examplecpp);
examplecppqd: $(examplecppqd);
examplef: $(examplef);

default: examplecpp ;

$(examplec): main.c $(COLESO_LIB)
	@echo  compiling $@ ...
	@$(CCX) $(COMPILER_OPTIONS) $(INCLUDES) main.c $(COLESO_LIB) -lstdc++ -lm -o $@

$(examplecpp): main.cpp $(COLESO_LIB)
	@echo  compiling $@ ...
	@$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) main.cpp $(COLESO_LIB) -o $@

$(examplecppqd): mainqd.cpp $(COLESO_LIBQD)
	@echo  compiling $@ ...
	@$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) -I.. mainqd.cpp $(COLESO_LIBQD) $(OUTDIR)/libqd.a -o $@

$(examplef): main.f $(COLESO_LIB)
	@echo  compiling $@ ...
	@$(FCX) main.f $(INCLUDES) $(COLESO_LIB) -J $(OUTDIR) -lstdc++ -o $(OUTDIR)/$@

clean:
	rm -f $(examplec) $(examplecpp) $(examplecppqd) $(examplef)
