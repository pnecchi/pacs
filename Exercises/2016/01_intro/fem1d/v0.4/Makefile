#-------------------------------------------------------
# Makefile template with src, obj, include directories
#
# Author: Pierpaolo Necchi
# Email: pierpaolo.necchi@gmail.com
#
# Changelog:
# 	2016-03-19 - First version
#-------------------------------------------------------

# Target
TARGET   = main

#Compiler and flags
CXX      = g++ 			# Alternative: clang++
CPPFLAGS = -std=c++11
CXXFLAGS = -g -Wall -O0

# Directories
BINDIR   = bin
OBJDIR   = obj
SRCDIR   = src
HDRDIR   = include

# Files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
HEADERS = $(wildcard $(HDRDIR)/*.h)

# Derive objects from sources
OBJECTS := $(patsubst $(SRCDIR)/%,%,$(SOURCES))
OBJECTS := $(patsubst %.cpp,%.o,$(OBJECTS))
OBJECTS := $(addprefix $(OBJDIR)/,$(OBJECTS))

# Remove command
RM = rm -f

# Rules and Actions
.PHONY : all clean

all : $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET) : $(OBJECTS)
	$(CXX) $(LFLAGS) $(OBJECTS) -o $@
	@echo ">> Linked "$<" successfully"

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ -I$(HDRDIR)
	@echo ">> Compiled "$<" successfully"

clean :
	$(RM) $(OBJECTS)





