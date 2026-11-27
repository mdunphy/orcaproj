CC := gcc

# Common flags shared by all compilations
CFLAGS := -O2 -g -fPIC -march=native -funroll-loops -fno-plt
CFLAGS += -Wall -Wpedantic -Wextra -Walloca -Wvla
#CFLAGS += -Wenum-conversion # if gcc >= 14
#CFLAGS += -Ofast
#CFLAGS += -fno-inline
CFLAGS += -fno-omit-frame-pointer
#CFLAGS += -fopenmp
#CFLAGS += -fopt-info-vec -fdump-tree-vect-details

CPPFLAGS := -Iinclude
LIBS := -lm

# C99 for the library
CFLAGS_LIB := $(CFLAGS) -std=gnu99

# GNU99 + POSIX timers for the benchmark tools
CFLAGS_BIN := $(CFLAGS) -std=gnu99 -D_POSIX_C_SOURCE=199309L

# For testing linking to C++
CXX := g++
CXXFLAGS := $(CFLAGS) -std=c++17

# Project layout
SRC_DIR := src
INC_DIR := include
OBJ_DIR := obj
BIN_DIR := bin
LIB_DIR := lib

# Coefficients path
ORCA05_DATA_PATH ?= $(CURDIR)/data
CFLAGS_LIB += -DORCA05_DATA_PATH=\"$(ORCA05_DATA_PATH)\"

# C sources
SRC_LIB := \
    chebcheb.c chebcheb_data.c chebcheb_tables.c \
    pchip.c pchip_data.c pchip_tables.c \
    blended.c forward.c inverse.c theta.c \
    functions.c data.c orca_init.c solver.c xy_to_ej.c

SRC_BIN := bench.c
SRC_BIN_CPP := bench_cpp.cpp

# Derive object lists
OBJ_LIB := $(addprefix $(OBJ_DIR)/,$(SRC_LIB:.c=.o))
OBJ_BIN := $(addprefix $(OBJ_DIR)/,$(SRC_BIN:.c=.o))
OBJ_BIN_CPP := $(addprefix $(OBJ_DIR)/,$(SRC_BIN_CPP:.cpp=.o))

# Output artifacts
SHLIB := $(LIB_DIR)/liborca.so
BINARIES := $(addprefix $(BIN_DIR)/,$(SRC_BIN:.c=))
BINARIES_CPP := $(addprefix $(BIN_DIR)/,$(SRC_BIN_CPP:.cpp=))

# Dependency generation
DEPFLAGS := -MMD -MP
DEPS := $(OBJ_LIB:.o=.d) $(OBJ_BIN:.o=.d)
DEPS += $(OBJ_BIN_CPP:.o=.d)

# Default target
.PHONY: all
all: $(SHLIB) $(BINARIES) $(BINARIES_CPP)

# Ensure directories exist
$(OBJ_DIR) $(BIN_DIR) $(LIB_DIR):
	@mkdir -p $@

# Compile library objects
$(OBJ_LIB): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS_LIB) $(DEPFLAGS) -c $< -o $@

# Compile benchmark objects
$(OBJ_BIN): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS_BIN) $(DEPFLAGS) -c $< -o $@

# Compile C++ objects
$(OBJ_BIN_CPP): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

# Compile C++ test
$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SRC_DIR)/%.cpp $(SHLIB) | $(BIN_DIR)
	$(CXX) $< \
	      -L$(LIB_DIR) -Wl,-rpath,'$$ORIGIN/../lib' \
	      $(LDFLAGS) $(LIBS) -lorca -o $@

# Link shared library
$(SHLIB): $(OBJ_LIB) | $(LIB_DIR)
	$(CC) -shared -g $(CFLAGS_LIB) $(OBJ_LIB) $(LDFLAGS) $(LIBS) -o $@

# Link benchmark executables
$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(SRC_DIR)/%.c $(SHLIB) | $(BIN_DIR)
	$(CC) $(CFLAGS_BIN) $< \
	      -L$(LIB_DIR) -Wl,-rpath,'$$ORIGIN/../lib' \
	      $(LDFLAGS) $(LIBS) -lorca -o $@

.PHONY: clean distclean run
clean:
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d

distclean: clean
	$(RM) $(SHLIB) $(BINARIES) $(BINARIES_CPP)
	-@rmdir $(OBJ_DIR) $(BIN_DIR) $(LIB_DIR) 2>/dev/null || true

run: $(BIN_DIR)/bench
	$(BIN_DIR)/bench

# Include auto-generated header dependencies
-include $(DEPS)
