# Compilador y flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
LIBS = -lboost_system -lpthread

# Directorios
SRC_DIR = .
OBJ_DIR = obj
BIN_DIR = bin

# Archivos
TARGET = $(BIN_DIR)/main
TARGET_K = $(BIN_DIR)/main_k
SOURCES = main.cpp
SOURCES_K = main_k.cpp
OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)
OBJECTS_K = $(SOURCES_K:%.cpp=$(OBJ_DIR)/%.o)

# Regla por defecto
all: directories $(TARGET) $(TARGET_K)

# Crear directorios necesarios
directories:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR)

# Compilar el ejecutable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Compilar ejecutable para k queries
$(TARGET_K): $(OBJECTS_K)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Compilar archivos objeto
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ejecutar el programa
run: all
	./$(TARGET)

# Ejecutar con k queries
run-k: all
	./$(TARGET_K)

# Limpiar archivos compilados
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Limpiar y recompilar
rebuild: clean all

.PHONY: all directories run clean rebuild
