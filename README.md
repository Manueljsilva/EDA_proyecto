# DB-LSH Demo (L=1, K=2, 2D) con R*-Tree

Pequeña demo inspirada en **DB-LSH** (Dynamic Bucket LSH): indexa proyecciones LSH de puntos 2D en un R*-Tree y responde consultas aproximadas mediante window queries sobre el espacio proyectado. La verificación final se hace con distancia real en el espacio original.

## 📋 Descripción

Este proyecto implementa una versión didáctica de DB-LSH (Database-friendly Locality-Sensitive Hashing) que combina:

- **Proyecciones LSH**: Transforma puntos 2D usando funciones hash lineales
- **R*-Tree indexing**: Almacena las proyecciones para búsquedas espaciales eficientes
- **Window queries dinámicas**: Búsqueda aproximada con expansión adaptativa
- **Algoritmos del paper**: Implementa `(r,c)-NN` y `c-ANN` queries

### Parámetros actuales

- **L = 1**: Una tabla hash (versión simplificada)
- **K = 2**: Dos funciones hash (proyecciones 2D)
- **C = 1.5**: Factor de aproximación
- **t = 1**: Parámetro de tolerancia
- **w₀ = 9.0**: Ancho base de ventana (4·C²)

### Funciones hash

```
h₁(x,y) = 0.6·x + 0.8·y
h₂(x,y) = 0.3·x - 0.9·y
```

## 🚀 Cómo ejecutar

### Requisitos

- **C++17** o superior
- **Boost.Geometry** (para R*-tree)
- **Make**

### Instalación de dependencias

```bash
# Ubuntu/Debian
sudo apt-get install libboost-all-dev

# Arch Linux
sudo pacman -S boost

# macOS (con Homebrew)
brew install boost
```

### Compilación y ejecución

```bash
# Compilar
make

# Ejecutar
./bin/main

# Limpiar y recompilar
make clean
make
```

### Comandos del Makefile

```bash
make          # Compilar el proyecto
make clean    # Limpiar archivos objeto y binarios
make rebuild  # Limpiar y recompilar
make run      # Compilar y ejecutar
```

## 📂 Estructura del proyecto

```
V_1/
├── main.cpp           # Implementación completa de DB-LSH
├── R_star.h           # Interfaz del R*-tree
├── R_star.cpp         # Implementación del R*-tree
├── Makefile           # Sistema de compilación
├── README.md          # Este archivo
├── bin/               # Ejecutables (generado)
├── obj/               # Archivos objeto (generado)
└── Dataset/
    └── Por verse 
```

## 🔍 Ejemplo de uso

```cpp
// Dataset de ejemplo
vector<tuple<double, double>> datos = {
    {1.0, 1.0}, {2.0, 2.0}, {4.0, 2.0},
    {5.0, 5.0}, {7.0, 8.0}
};

// Crear índice DB-LSH
DBfsh indice(2, 1, 2, 1.5, 1);
indice.insertar(datos);

// Consulta c-ANN
tuple<double, double> query = {6.0, 6.0};
tuple<double,double> vecino = indice.C_ANN(query, 1.5);
```

### Salida esperada

```
=====================================
DB-LSH con (r,c)-NN y c-ANN Query
=====================================
Hashes generados:
Original: (1, 1) -> Hash: (1.4, -0.6)
Original: (5, 5) -> Hash: (7, -3)
...

c-ANN Query
Query: (6, 6)
c = 1.5

--- (r,c)-NN Query ---
Ventana W(G(q), w_0·r = 9): [3.9, 12.9] x [-8.1, 0.9]
✓ Condición ||q,o|| ≤ cr cumplida

RESULTADO FINAL:
  Vecino más cercano (c-aproximado): (5, 5)
  Distancia real: 1.41421
```

## 🧪 Algoritmos implementados

### Algorithm 1: (r,c)-NN Query

Busca un vecino aproximado dentro de un radio `r` con factor de aproximación `c`:

- Calcula hash del query: `G(q)`
- Define ventana: `W(G(q), w₀·r)`
- Busca candidatos en el R*-tree
- Retorna si `||q,o|| ≤ c·r` o si `cnt = 2tL + 1`

### Algorithm 2: c-ANN Query

Encuentra un vecino c-aproximado con expansión dinámica:

- Comienza con `r = 1`
- Llama a `(r,c)-NN`
- Si no encuentra, expande `r ← c·r`
- Repite hasta encontrar un vecino

## 🎯 Mejoras futuras

### Corto plazo

- ⬜ **Ventanas dinámicas**: `w = w₀ · rᵢ` con expansión `rᵢ₊₁ = c · rᵢ` (paper-like)
- ⬜ **Deduplicación de candidatos**: `unordered_set<int>` entre iteraciones
- ⬜ **a aleatorios**: `𝒩(0,1)` normalizados, con semilla fija para reproducibilidad
- ⬜ **Separar datos**: Índice guarda solo proyecciones + id; vector aparte con originales
- ⬜ **CSV (fstream)**: Lector robusto con `reserve()`, `from_chars`, etc. y validación
- ⬜ **Métricas**: Imprimir `recall@k`, `overall ratio`, #candidatos, tiempo por query

### Medio plazo

- ⬜ **k-ANN**: Con `priority_queue` (max-heap tamaño k)
- ⬜ **Multi-L**: Soportar varios grupos de proyección (sube recall)
- ⬜ **Parámetros por CLI**: `--w0 --c --r0 --L --K --seed`
- ⬜ **Tests unitarios**: Básicos (proyección, ventana, verificación)

### Largo plazo

- ⬜ **Multi-probe**: Ventanas adyacentes priorizadas
- ⬜ **Alto-dimensional**: (d>2), soporte genérico K-D en Boost
- ⬜ **Persistencia**: Del índice y batch de queries
- ⬜ **Benchmarks**: Con datasets públicos (SIFT, GloVe, Deep, etc.)

## 📊 Características técnicas

| Característica | Valor |
|----------------|-------|
| **Dimensiones** | 2D (fijo) |
| **Tablas hash** | L = 1 |
| **Funciones hash** | K = 2 |
| **Estructura espacial** | R*-tree (Boost) |
| **Parámetro R*-tree** | 16 elementos/nodo |
| **Tipo de consulta** | c-ANN aproximada |
| **Verificación** | Distancia euclidiana real |

## 🔧 Detalles de implementación

### Clase DBfsh

```cpp
class DBfsh {
private:
    int D;              // Dimensiones (2)
    int L;              // Tablas hash (1)
    int K;              // Funciones hash (2)
    double C;           // Factor aproximación (1.5)
    double w0;          // Ancho ventana base (9.0)
    int t;              // Parámetro tolerancia (1)
    RStarTreeIndex indice;  // R*-tree de Boost
    
public:
    void insertar(vector<tuple<double,double>> datos);
    tuple<double,double> RC_NN(tuple<double,double> q, double r, double c);
    tuple<double,double> C_ANN(tuple<double,double> q, double c);
};
```

### Mapeo hash ↔ original

El sistema mantiene un mapeo bidireccional:
- **Inserción**: `original → hash` (proyección)
- **Búsqueda**: `hash → original` (recuperación)

Esto permite trabajar en el espacio proyectado pero retornar puntos originales.

## 📖 Referencias

Este proyecto está inspirado en técnicas de **Locality-Sensitive Hashing** para búsqueda aproximada de vecinos más cercanos (ANN) en espacios de alta dimensión.

### Conceptos clave

- **LSH**: Proyecta puntos similares a buckets similares
- **DB-LSH**: Variante sin buckets fijos, usa window queries
- **R*-tree**: Índice espacial para consultas rectangulares eficientes
- **c-ANN**: Encuentra vecinos a distancia ≤ c·dist(q, NN)

## 🤝 Contribuciones

Este es un proyecto educativo. Sugerencias y mejoras son bienvenidas.

## 📄 Licencia

Proyecto académico para Estructuras de Datos Avanzadas.

---

**Autor**: Manuel J. Simpson  
**Fecha**: Noviembre 2025  
**Tecnologías**: C++17, Boost.Geometry, R*-tree, LSH
