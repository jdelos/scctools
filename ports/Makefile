# Compiler and Flags
EMCC ?= emcc
CXXFLAGS ?= -std=c++17 -O0 -Wall -Wextra -g3 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1

# Include and Library Paths. Keep SymEngine headers and library from same install.
# SymEngine 0.13 headers and library must match; override for native/WASM installs.
SYMENGINE_INCLUDE ?= /usr/local/include
SYMENGINE_LIB ?= /usr/local/lib
INCLUDES = -Iinclude \
           -I$(SYMENGINE_INCLUDE) \
           -I./libs/mpfr-4.2.1/install/include \
           -I./libs/gmp-6.2.1/install/include

LDFLAGS = -L./libs/symengine/build_wasm/symengine/ \
          -L./libs/mpfr-4.2.1/install/lib \
          -L./libs/gmp-6.2.1/install/lib \
          -lsymengine -lgmp -lmpfr

# Output files
TARGET = test_dickson_matrix
WASM_DIR ?= wasm
TEST_SRC_DICKSON_MATRIX ?= test/test_dickson_matrix.cpp
TEST_TARGET_DICKSON_MATRIX = $(WASM_DIR)/test_dickson_matrix.js

EXPORT_FUNCS = -s EXPORTED_FUNCTIONS="['_main']"
EXPORT_RUNTIME = -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']"
ENV_SETTINGS = -s ENVIRONMENT='web'

# Source and Object Files
SRC_DIR = src
OBJ_DIR = obj
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp) main.cpp
OBJ_FILES = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES)))

# Rules
all: $(TARGET)

native-test:
	g++ -std=c++17 -Wall -Wextra -Iinclude -I$(SYMENGINE_INCLUDE) test/test_dickson_arch.cpp src/dickson_arch.cpp src/dickson_matrix.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_dickson_arch
	/tmp/test_dickson_arch

native-phase-test:
	g++ -std=c++17 -Wall -Wextra -Iinclude -I$(SYMENGINE_INCLUDE) test/test_phase.cpp src/dickson_hybrid_topology.cpp src/dickson_arch.cpp src/dickson_matrix.cpp src/phase_conv.cpp src/scc_phase.cpp src/graph_primitives.cpp src/utilities.cpp src/solve_charge_vectors.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_phase
	/tmp/test_phase

native-charge-test:
	g++ -std=c++17 -Wall -Wextra -Iinclude -I$(SYMENGINE_INCLUDE) test/test_solve_charge_vectors.cpp src/solve_charge_vectors.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_solve_charge_vectors
	/tmp/test_solve_charge_vectors

native-graph-primitives-test:
	g++ -std=c++17 -Wall -Wextra -Iinclude -I$(SYMENGINE_INCLUDE) test/test_graph_primitives.cpp src/graph_primitives.cpp -L$(SYMENGINE_LIB) -lsymengine -lgmp -o /tmp/test_graph_primitives
	/tmp/test_graph_primitives

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(EMCC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/main.o: main.cpp | $(OBJ_DIR)
	$(EMCC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link and generate WebAssembly output
$(TARGET): $(OBJ_FILES)
	$(EMCC) $(CXXFLAGS) $(INCLUDES) $(OBJ_FILES) $(LDFLAGS) \
		-o $(TARGET) $(EXPORT_FUNCS) $(EXPORT_RUNTIME) $(ENV_SETTINGS)

# Individual Function Test - dickson_matrix
test_dickson_matrix: $(TEST_TARGET_DICKSON_MATRIX)

$(TEST_TARGET_DICKSON_MATRIX): $(TEST_SRC_DICKSON_MATRIX) | $(OBJ_DIR)
	$(EMCC) $(CXXFLAGS) $(INCLUDES) $(TEST_SRC_DICKSON_MATRIX) -o $(TEST_TARGET_DICKSON_MATRIX) $(EXPORT_FUNCS) $(EXPORT_RUNTIME)
	@echo "Test build completed: $(TEST_TARGET_DICKSON_MATRIX)"
	@echo "Run by opening wasm/test_runner.html in the browser."

# Create object directory if not exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean build files
clean:
	rm -rf $(OBJ_DIR) $(TARGET)
