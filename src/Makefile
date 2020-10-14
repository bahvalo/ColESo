CXX?=mpic++
COMPILER_OPTIONS = -O3
INCLUDES = -I. -I../
LIBS ?=
QDSUFF ?=
OUTDIR ?= ../bin

#make CFLAGS=-DEXTRAPRECISION_COLESO QDSUFF=qd

TARGET = $(addprefix $(OUTDIR)/, coleso$(QDSUFF).a)

SOURCES = $(wildcard es_*.cpp) io_parser.cpp
BINOBJ = $(addprefix $(OUTDIR)/$(QDSUFF), $(patsubst %.cpp, %.o, $(SOURCES)))

default: $(TARGET) ;

$(OUTDIR)/$(QDSUFF)%.o : %.cpp
	$(CXX) $(COMPILER_OPTIONS) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(BINOBJ): | $(OUTDIR)

$(TARGET): $(BINOBJ) 
	ar rv $@ $(BINOBJ)
	ranlib $@

clean:
	rm -f $(BINOBJ) $(TARGET)
