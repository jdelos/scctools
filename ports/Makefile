# Native production build. Override SYMENGINE_INCLUDE and SYMENGINE_LIB for other installs.
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O0 -Wall -Wextra -g3
SYMENGINE_INCLUDE ?= /usr/local/include
SYMENGINE_LIB ?= /usr/local/lib
INCLUDES = -Iinclude -I$(SYMENGINE_INCLUDE)
LDFLAGS = -L$(SYMENGINE_LIB) -lsymengine -lgmp -lmpfr

TARGET = scctools
SRC_DIR = src
OBJ_DIR = obj
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp) main.cpp
OBJ_FILES = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES)))
TEST_SRC_DICKSON_MATRIX ?= test/test_dickson_matrix.cpp

.PHONY: all native-test native-phase-test native-charge-test native-graph-primitives-test native-matrix-test native-boundary-test wasm-test wasm-build test_dickson_matrix clean FORCE

all: $(TARGET)

wasm-build:
	./wasm/build.sh

wasm-test: wasm-build
	node test/test_wasm.js

native-test:
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/test_dickson_arch.cpp src/dickson_arch.cpp src/dickson_matrix.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_dickson_arch
	/tmp/test_dickson_arch

native-phase-test:
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/test_phase.cpp src/dickson_hybrid_topology.cpp src/dickson_arch.cpp src/dickson_matrix.cpp src/phase_conv.cpp src/scc_phase.cpp src/graph_primitives.cpp src/utilities.cpp src/solve_charge_vectors.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_phase
	/tmp/test_phase

native-charge-test:
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/test_solve_charge_vectors.cpp src/solve_charge_vectors.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_solve_charge_vectors
	/tmp/test_solve_charge_vectors

native-graph-primitives-test:
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/test_graph_primitives.cpp src/graph_primitives.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_graph_primitives
	/tmp/test_graph_primitives

native-matrix-test:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TEST_SRC_DICKSON_MATRIX) src/dickson_matrix.cpp -o /tmp/test_dickson_matrix
	/tmp/test_dickson_matrix

native-boundary-test:
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/test_native_boundary.cpp src/native_boundary.cpp src/dickson_hybrid_topology.cpp src/dickson_arch.cpp src/dickson_matrix.cpp src/phase_conv.cpp src/scc_phase.cpp src/graph_primitives.cpp src/utilities.cpp src/solve_charge_vectors.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_native_boundary
	/tmp/test_native_boundary

test_dickson_matrix: native-matrix-test

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp FORCE | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/main.o: main.cpp FORCE | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ_FILES) $(LDFLAGS) -o $@

FORCE:

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
