# Makefile for SEED 4.0.0
# clangini 2016
#SHELL = /bin/bash
# Folders:
#SRC_DIR = .
BUILD_DIR = build
BUILD_DIR_MPI = build/mpi
TARGET_EXEC = seed_4
TARGET_EXEC_MPI = seed_4_mpi

# Flags:
CXX = g++
CXXFLAGS = -std=c++11 -O3 -pedantic -Wall -march=native
LDFLAGS = -lm -O3 -static

INC_FLAGS = -I. -I./lib/

# For seed MPI:
MPICXX = mpic++
MPICXXFLAGS = -std=c++11 -O3 -pedantic -Wall -march=native
MPILDFLAGS = -lm -O3

include Makefile.local

DEP_FLAGS = $(INC_FLAGS) -MMD -MP

PRESRCLIST = align_frag assire checkresn cofren cogrep conubu convert2pp count_heavy_atom extoutnam featres fialhy findfrsym fracdo geomcent geomfu hybstat listchres main makbsatlist makspv nrutil otfunc psspee psspma reacdo reduc_polvectre refrfi_coo_mol2 refrfi_mol2_new reinfi relaic relaic_en rerefi_mol2 rosefr seedfr_ap seedfr simila solv_frag solv_frag_fast solv_lookup solv_util sort_all sort sqdisfrre_ps vwfren vwgrep write_charmm write_chm_clus write_mol2
SRCLIST = $(PRESRCLIST:%=%.cpp)
SRCS = $(SRCLIST)
OBJS = $(SRCLIST:%.cpp=$(BUILD_DIR)/%.o)
OBJS_MPI = $(SRCLIST:%.cpp=$(BUILD_DIR_MPI)/%.o)
DEPS = $(OBJS:%.o=%.d)
DEPS_MPI = $(OBJS_MPI:%.o=%.d)

# Commands:
RM = rm
MKDIR_P = mkdir -p
RMDIR = rmdir
ECHO = @echo

# Targets:
.PHONY: seed
seed: $(TARGET_EXEC)

.PHONY: seed_4_mpi
seed_mpi: $(TARGET_EXEC_MPI)

$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

$(TARGET_EXEC_MPI): $(OBJS_MPI)
	$(MPICXX) $(OBJS_MPI) $(MPILDFLAGS) -o $@

# Using order-only prerequisite to build the directory
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) -c $(CXXFLAGS) $(DEP_FLAGS) $< -o $@

$(BUILD_DIR_MPI)/%.o: %.cpp | $(BUILD_DIR_MPI)
	$(MPICXX) -c $(MPICXXFLAGS) -DENABLE_MPI $(DEP_FLAGS) $< -o $@

$(BUILD_DIR):
	$(MKDIR_P) $@

$(BUILD_DIR_MPI): $(BUILD_DIR)
	$(MKDIR_P) $@

-include $(DEPS)
-include $(DEPS_MPI)

.PHONY: print
print:
	$(ECHO) $(PRESRCLIST)
	$(ECHO) $(SRCLIST)
	$(ECHO) $(SRCS)
	$(ECHO) $(OBJS)
	$(ECHO) $(DEPS)

.PHONY: objclean
objclean:
	for i in $(OBJS_MPI); do (if [ -e $${i} ]; then $(RM) $${i}; fi;) done
	for i in $(DEPS_MPI); do (if [ -e $${i} ]; then $(RM) $${i}; fi;) done
	for i in $(OBJS); do (if [ -e $${i} ]; then $(RM) $${i}; fi;) done
	for i in $(DEPS); do (if [ -e $${i} ]; then $(RM) $${i}; fi;) done
	if [ -d $(BUILD_DIR_MPI) ]; then $(RMDIR) $(BUILD_DIR_MPI); fi
	if [ -d $(BUILD_DIR) ]; then $(RMDIR) $(BUILD_DIR); fi

.PHONY: clean
clean: objclean
	if [ -e $(TARGET_EXEC) ]; then $(RM) $(TARGET_EXEC); fi
	if [ -e $(TARGET_EXEC_MPI) ]; then $(RM) $(TARGET_EXEC_MPI); fi
