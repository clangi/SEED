# Makefile for SEED 4.0.0
# clangini 2016
#SHELL = /bin/bash
# Folders:
#SRC_DIR = .
BUILD_DIR = build
TARGET_EXEC = seed_4_linesearch

# Flags:
CXXFLAGS = -std=c++11 -O3 -pedantic -Wall -march=native
CXX = g++
LDFLAGS = -lm -O3 -static
INC_FLAGS = -I.

include Makefile.local

DEP_FLAGS = $(INC_FLAGS) -MMD -MP

PRESRCLIST = align_frag assire checkresn cofren cogrep conubu convert2pp count_heavy_atom extoutnam featres fialhy findfrsym fracdo geomcent geomfu hybstat listchres main makbsatlist makspv nrutil otfunc pose psspee psspma reacdo reduc_polvectre refrfi_coo_mol2 refrfi_mol2_new reinfi relaic relaic_en rerefi_mol2 rosefr seedfr_ap seedfr simila solv_frag solv_frag_fast solv_lookup solv_util sort_all sort sqdisfrre_ps vwfren vwgrep write_charmm write_chm_clus write_mol2 
SRCLIST = $(PRESRCLIST:%=%.c)
#SRCS = $(SRCLIST:%=$(SRC_DIR)/%)
SRCS = $(SRCLIST)
OBJS = $(SRCLIST:%.c=$(BUILD_DIR)/%.o)
DEPS = $(OBJS:%.o=%.d)

# Commands:
RM = rm
MKDIR_P = mkdir -p
RMDIR = rmdir
ECHO = @echo

# Targets:
seed: $(TARGET_EXEC)

$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	$(MKDIR_P) $(BUILD_DIR)
	$(CXX) -c $(CXXFLAGS) $(DEP_FLAGS) $< -o $@

-include $(DEPS)

.PHONY: print
print:
	$(ECHO) $(PRESRCLIST)
	$(ECHO) $(SRCLIST)
	$(ECHO) $(SRCS)
	$(ECHO) $(OBJS)
	$(ECHO) $(DEPS)

.PHONY: clean
clean:
	if [ -e $(TARGET_EXEC) ]; then $(RM) $(TARGET_EXEC); fi
	for i in $(OBJS); do (if [ -e $${i} ]; then $(RM) $${i}; fi;) done
	for i in $(DEPS); do (if [ -e $${i} ]; then $(RM) $${i}; fi;) done
	if [ -d $(BUILD_DIR) ]; then $(RMDIR) $(BUILD_DIR); fi
