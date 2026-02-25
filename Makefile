ROOTCONFIG = root-config
ROOTCFLAGS = $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS   = $(shell $(ROOTCONFIG) --libs)

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra -I./include $(ROOTCFLAGS)

SRCDIR = src
OBJDIR = obj
BINDIR = .

LIB_SOURCES = $(wildcard $(SRCDIR)/*.cxx)
LIB_OBJECTS = $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(LIB_SOURCES))

TARGET = $(BINDIR)/bwgen

.PHONY: all clean

all: $(OBJDIR) $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(TARGET): $(LIB_OBJECTS) apps/bwgen.cpp
	$(CXX) $(CXXFLAGS) -o $@ apps/bwgen.cpp $(LIB_OBJECTS) $(ROOTLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) $(TARGET)
