# Compiler and Flags
EMCC = emcc
CXXFLAGS = -std=c++17 -O0 -Wall -Wextra -g3 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1

# Include and Library Paths
INCLUDES = -Iinclude \
           -I/usr/local/include \
           -I./libs/mpfr-4.2.1/install/include \
           -I./libs/gmp-6.2.1/install/include

LDFLAGS = -L./libs/symengine/build_wasm/symengine/ \
          -L./libs/mpfr-4.2.1/install/lib \
          -L./libs/gmp-6.2.1/install/lib \
          -lsymengine -lgmp -lmpfr

# Output files
TARGET = test_dickson_matrix
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
