######################################################################
# Makefile
######################################################################

######################################################################
# where to find the .h files
VPATH =	utilities:smt

# list ot the source files (with subdirectory, without .cc extension)
FILES = utilities/utilities \
	smt/configuration smt/paths \
	smt/cc-dubins smt/cc-rs \
	main/test

# name of the bin file
BIN   = ccsteer

# name of the doxygen configuration file
DOXYFILE = Doxyfile

######################################################################

DIRS		:= $(subst :, ,$(VPATH))
INCL_DIRS	:= $(DIRS:%=-I%)

CC_FILES	:= $(FILES:%=%.cc)
O_FILES		:= $(FILES:%=%.o)

CXXFLAGS = -g -Wall $(INCL_DIRS)
LDLIBS = 

# dependencies file
DEPFILE = .dependencies
FOODEPFILE := $(shell touch $(DEPFILE))

######################################################################
.PHONY: all
all: depend bin

######################################################################
# create the executable file
.PHONY: bin
bin: $(O_FILES)
	$(CXX) $(CXXFLAGS) -o $(BIN) $^ $(LDLIBS) 

######################################################################
# remove the object files and the executable file
.PHONY: clean
clean:
	${RM} $(O_FILES) ${BIN}
	find . \( -name core -o -name "*~" -o -name "#*#" -o -name ".*~" \) -exec ${RM} \{} \;

######################################################################
# update the doxygen documentation
.PHONY: doxygen
doxygen:
	doxygen ${DOXYFILE}

######################################################################
# update the dependencies file
.PHONY: depend
depend:
	echo --- updating $(DEPFILE)
	makedepend 	$(INCL_DIRS) -Y \
			$(CC_FILES) \
			-f $(DEPFILE)

include $(DEPFILE)
# DO NOT DELETE
