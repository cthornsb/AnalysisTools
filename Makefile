CC = g++

CFLAGS = -g -fPIC -Wall -O3 `root-config --cflags --glibs` -Iinclude
LDLIBS = -lm -lstdc++ -lgsl -lgslcblas -lgfortran `root-config --libs`
ROOT_INC = `root-config --incdir`

TOP_LEVEL = $(shell pwd)
DICT_DIR = $(TOP_LEVEL)/dict
INCLUDE_DIR = $(TOP_LEVEL)/include
SOURCE_DIR = $(TOP_LEVEL)/src

OBJ_DIR = $(TOP_LEVEL)/obj
DICT_OBJ_DIR = $(DICT_DIR)/obj

TOOLS = TimeAlign TimeAlignCf CheckAlign PulseViewer PulseAnalyzer Stitcher Gater Overlay

# ROOT dictionary stuff
DICT_SOURCE = RootDict
STRUCT_FILE = Structures

ROOTOBJ = $(DICT_OBJ_DIR)/$(DICT_SOURCE).o 
ROOTOBJ += $(OBJ_DIR)/$(STRUCT_FILE).o
SFLAGS = $(addprefix -l,$(DICT_SOURCE))

#####################################################################

all: directory $(DICT_OBJ_DIR)/$(DICT_SOURCE).so $(TOOLS)
#	Create everything we need to make tools

dictionary: $(DICT_OBJ_DIR)/$(DICT_SOURCE).so
#	Create root dictionary objects

.PHONY: clean tidy directory

.SECONDARY: $(DICT_DIR)/$(DICT_SOURCE).cpp $(ROOTOBJ)
#	Want to keep the source files created by rootcint after compilation
#	as well as keeping the object file made from those source files

#####################################################################

directory: $(OBJ_DIR) $(DICT_OBJ_DIR)

$(OBJ_DIR):
#	Make the object file directory
	mkdir $(OBJ_DIR)

$(DICT_OBJ_DIR):
#	Make root dictionary object file directory
	mkdir $(DICT_OBJ_DIR)

########################################################################

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
#	Compile C++ source files
	$(CC) -c $(CFLAGS) $< -o $@

#####################################################################

$(DICT_OBJ_DIR)/%.o: $(DICT_DIR)/%.cpp
#	Compile rootcint source files
	$(CC) -c $(CFLAGS) $< -o $@

$(DICT_OBJ_DIR)/%.so: $(DICT_OBJ_DIR) $(OBJ_DIR)/Structures.o $(DICT_OBJ_DIR)/$(DICT_SOURCE).o
#	Generate the root shared library (.so) for the dictionary
	$(CC) -g -shared -Wl,-soname,lib$(DICT_SOURCE).so -o $(DICT_OBJ_DIR)/lib$(DICT_SOURCE).so $(OBJ_DIR)/Structures.o $(DICT_OBJ_DIR)/$(DICT_SOURCE).o -lc

$(DICT_DIR)/%.cpp: $(INCLUDE_DIR)/$(STRUCT_FILE).h $(DICT_DIR)/LinkDef.h
#	Generate the dictionary source files using rootcint
	@cd $(DICT_DIR); rootcint -f $@ -c $(INCLUDE_DIR)/$(STRUCT_FILE).h $(DICT_DIR)/LinkDef.h

#####################################################################

TimeAlign: $(SOURCE_DIR)/TimeAlign.cpp
#	Compile TimeAlign tool
	$(CC) $(CFLAGS) $< -o $@

TimeAlignCf: $(SOURCE_DIR)/TimeAlignCf.cpp
#	Compile TimeAlignCf tool
	$(CC) $(CFLAGS) $< -o $@
	
CheckAlign: $(SOURCE_DIR)/CheckAlign.cpp
#	Compile CheckAlign tool
	$(CC) $(CFLAGS) $< -o $@

PulseViewer: $(SOURCE_DIR)/PulseViewer.cpp
#	Compile PulseViewer tool
	$(CC) $(CFLAGS) $< -o $@

PulseAnalyzer: $(SOURCE_DIR)/PulseAnalyzer.cpp
#	Compile PulseAnalyzer tool
	$(CC) $(CFLAGS) $< -o $@
	
Stitcher: $(SOURCE_DIR)/Stitcher.cpp
#	Compile Stitcher tool
	$(CC) $(CFLAGS) $< -o $@
	
Gater: $(DICT_OBJ_DIR)/$(DICT_SOURCE).so $(OBJ_DIR)/Gater.o
#	Compile Gater tool
	$(CC) $(OBJ_DIR)/$@.o $(ROOTOBJ) -L$(DICT_OBJ_DIR) $(SFLAGS) -o $@ $(LDLIBS)
	
Overlay: $(SOURCE_DIR)/Overlay.cpp
#	Compile Stitcher tool
	$(CC) $(CFLAGS) $< -o $@

#####################################################################

tidy: clean_obj

clean: clean_obj clean_dict

clean_obj:
	@echo "Cleaning up..."
	@rm -f $(OBJ_DIR)/*.o $(TOOLS)
	
clean_dict:
	@echo "Removing ROOT dictionaries..."
	@rm -f $(DICT_DIR)/$(DICT_SOURCE).cpp $(DICT_DIR)/$(DICT_SOURCE).h $(DICT_OBJ_DIR)/*.o  $(DICT_OBJ_DIR)/*.so
