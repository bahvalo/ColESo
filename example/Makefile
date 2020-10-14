CXX?=c++
CCX?=cc
FCX?=gfortran
COMPILER_OPTIONS = -O3
INCLUDES = -I../src
OUTDIR = ../bin
LIBS ?= $(OUTDIR)/coleso.a

default: examplecpp ;

examplec: $(OUTDIR)/examplec
$(OUTDIR)/examplec:
	@echo  compiling $@ ...
	@$(CCX) $(COMPILER_OPTIONS) $(INCLUDES) main.c $(LIBS) -lstdc++ -lm -o $@

examplecpp: $(OUTDIR)/examplecpp
$(OUTDIR)/examplecpp:
	@echo  compiling $@ ...
	@$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) main.cpp $(LIBS) -o $@

examplecppqd: $(OUTDIR)/examplecppqd
$(OUTDIR)/examplecppqd:
	@echo  compiling $@ ...
	@$(CXX) $(COMPILER_OPTIONS) $(INCLUDES) -I.. mainqd.cpp $(OUTDIR)/colesoqd.a $(OUTDIR)/libqd.a -o $@

examplef: $(OUTDIR)/examplef
$(OUTDIR)/examplef:
	@echo  compiling $@ ...
	@$(FCX) main.f $(INCLUDES) $(LIBS) -J $(OUTDIR) -lstdc++ -o $(OUTDIR)/$@

clean:
	rm -f $(OUTDIR)/testc $(OUTDIR)/testcpp $(OUTDIR)/testf
