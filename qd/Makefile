CXX?=c++

OPTFLAGS = -O2
INCLUDES = -I..
COPTIONS = -w
OUTDIR ?= ../bin

SRC = bits.cpp c_dd.cpp c_qd.cpp dd_const.cpp dd_real.cpp fpu.cpp qd_const.cpp qd_real.cpp util.cpp

TARGET := $(addprefix $(OUTDIR)/, libqd.a)
OBJS   := $(addsuffix .o,$(basename $(SRC)))
BINOBJ := $(addprefix $(OUTDIR)/, $(OBJS))


#CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
#LD = $(CC) -L. 

default: $(TARGET) ;

$(OUTDIR)/%.o : %.cpp
	$(CXX) $(COPTIONS) $(OPTFLAGS) $(INCLUDES) -c  $< -o $@

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(BINOBJ): | $(OUTDIR)

$(TARGET): $(BINOBJ) 
	ar rv $@ $(BINOBJ)
	ranlib $@

clean:
	rm -f $(BINOBJ) $(TARGET)
