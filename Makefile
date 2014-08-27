CC = g++

TOOL_SRC_DIR = $(shell pwd)/src

TIMEALIGN = TimeAlign
TIMEALIGNCF = TimeAlignCf
CHECKALIGN = CheckAlign
VIEWER = PulseViewer
ANALYZER = PulseAnalyzer
STITCHER = Stitcher
GATER = Gater

ALL_TOOLS = $(TIMEALIGN) $(TIMEALIGNCF) $(CHECKALIGN) $(VIEWER) $(ANALYZER) $(STITCHER) $(GATER)

all: $(ALL_TOOLS)
#	Make all tools

clean:
#	Clean up
	@rm -f $(ALL_TOOLS)

$(TIMEALIGN): $(TOOL_SRC_DIR)/$(TIMEALIGN).cpp
#	Make the TimeAlign tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(TIMEALIGN).cpp `root-config --cflags --glibs` -o $@
	
$(TIMEALIGNCF): $(TOOL_SRC_DIR)/$(TIMEALIGNCF).cpp
#	Make the TimeAlign tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(TIMEALIGNCF).cpp `root-config --cflags --glibs` -o $@

$(CHECKALIGN): $(TOOL_SRC_DIR)/$(CHECKALIGN).cpp
#	Make the CheckAlign tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(CHECKALIGN).cpp `root-config --cflags --glibs` -o $@
	
$(VIEWER): $(TOOL_SRC_DIR)/$(VIEWER).cpp
#	Make the PulseViewer tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(VIEWER).cpp `root-config --cflags --glibs` -o $@
	
$(ANALYZER): $(TOOL_SRC_DIR)/$(ANALYZER).cpp
#	Make the PulseAnalyzer tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(ANALYZER).cpp `root-config --cflags --glibs` -o $@

$(STITCHER): $(TOOL_SRC_DIR)/$(STITCHER).cpp
#	Make the Stitcher tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(STITCHER).cpp `root-config --cflags --glibs` -o $@
	
$(GATER): $(TOOL_SRC_DIR)/$(GATER).cpp
#	Make the Gater tool
	$(CC) -O2 -Wall $(TOOL_SRC_DIR)/$(GATER).cpp `root-config --cflags --glibs` -o $@
