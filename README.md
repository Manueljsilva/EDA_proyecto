# DB-LSH N-Dimensional (V3) con R*-Tree 10D

Implementación de **DB-LSH** (Database-friendly Locality-Sensitive Hashing) siguiendo el paper original, con extensiones para **N dimensiones** y proyección a **K=10 dimensiones**. Utiliza **window queries dinámicas** sobre un R*-tree 10D en lugar de buckets estáticos, permitiendo búsquedas aproximadas de vecinos más cercanos (c-ANN) eficientes en espacios de alta dimensionalidad.

## 📋 ¿Qué es DB-LSH?

**DB-LSH** es una variante de LSH diseñada para integrarse con estructuras de indexación espacial existentes (como R*-tree) en lugar de usar tablas hash tradicionales.

### Diferencias clave con LSH clásico:

| Característica | LSH clásico | DB-LSH (este proyecto) |
|----------------|-------------|------------------------|
| **Almacenamiento** | Buckets hash estáticos | R*-tree 10D dinámico |
| **Búsqueda** | Lookup directo en bucket | Window query 10D expansiva |
| **Colisiones** | Todos los puntos en bucket | Ventana `W(G(q), w₀·r)` en 10D |
| **Expansión** | Probar múltiples tablas L | Expandir radio `r ← c·r` |
| **Indexación** | Hash tables | Índice espacial R*-tree 10D |
| **Proyección** | Variable | N-D → 10D fijo |

### Ventajas de DB-LSH:

✅ **Sin buckets fijos**: No desperdicia memoria en buckets vacíos  
✅ **Consultas dinámicas**: Ajusta la ventana según necesidad (r expansivo)  
✅ **Integración DB**: Se integra con índices espaciales existentes  
✅ **Escalable**: Funciona bien con datasets grandes

## 🎯 Algoritmos del paper implementados

Este proyecto implementa fielmente los **Algorithm 1** y **Algorithm 2** del paper DB-LSH:

### **Algorithm 1: (r,c)-NN Query**

Busca un punto `o` tal que `||q,o|| ≤ c·r` (vecino c-aproximado dentro de radio r):

```
Input: q (query point), r (radio), c (approximation ratio), t (parámetro)
Output: punto o ó ∅

cnt ← 0
for i = 1 to L do:
    Compute G_i(q)  // proyección hash del query
    while a point o ∈ W(G_i(q), w₀·r) is found do:
        cnt ← cnt + 1
        if cnt = 2tL + 1 OR ||q,o|| ≤ c·r then:
            return o
return ∅
```

**Idea clave**: Usa **window query** `W(G(q), w₀·r)` en el espacio proyectado para encontrar candidatos.

### **Algorithm 2: c-ANN Query** 

Encuentra un vecino c-aproximado expandiendo el radio dinámicamente:

```
Input: q (query point), c (approximation ratio)
Output: punto o

r ← 1
while TRUE do:
    o ← call (r,c)-NN
    if o ≠ ∅ then:
        return o
    else:
        r ← c·r  // expandir radio
```

**Idea clave**: Si no encuentra vecinos con radio `r`, **expande a `r ← c·r`** y repite.

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
cd V_2

# Compilar
make

# Ejecutar
./bin/main

# O compilar y ejecutar directamente
make run

# Limpiar y recompilar
make clean && make
```

## 📂 Estructura del proyecto

```
V_3/
├── main.cpp           # Implementación K=10 con referencias al paper
├── R_star.h           # Interfaz del R*-tree 10D
├── R_star.cpp         # Implementación del R*-tree 10D
├── Makefile           # Sistema de compilación
├── README.md          # Este archivo
├── bin/               # Ejecutables (generado)
└── obj/               # Archivos objeto (generado)
```

## 🔍 Ejemplos de uso

### Ejemplo 1: 20D → 10D (reducción real)

```cpp
vector<vector<double>> datos_20d;
for(int i = 0; i < 8; i++) {
    vector<double> punto(20);
    for(int j = 0; j < 20; j++) {
        punto[j] = (i + 1) * 0.5 + j * 0.1;
    }
    datos_20d.push_back(punto);
}

DBfsh indice_20d(20, 1, 1.5, 1, 42);  // dim=20, L=1, C=1.5, t=1
indice_20d.insertar(datos_20d);

vector<double> query_20d(20);
for(int j = 0; j < 20; j++) {
    query_20d[j] = 1.5 + j * 0.1;
}
vector<double> vecino = indice_20d.C_ANN(query_20d, 1.5);
```

### Ejemplo 2: 50D → 10D

```cpp
vector<vector<double>> datos_50d;
mt19937 gen(123);
normal_distribution<double> dist(0.0, 1.0);

for(int i = 0; i < 12; i++) {
    vector<double> punto(50);
    for(int j = 0; j < 50; j++) {
        punto[j] = dist(gen) + i * 0.2;
    }
    datos_50d.push_back(punto);
}

DBfsh indice_50d(50, 1, 2.0, 1, 456);
indice_50d.insertar(datos_50d);

vector<double> query_50d = datos_50d[5];
vector<double> vecino = indice_50d.C_ANN(query_50d, 2.0);
```

### Ejemplo 3: 128D → 10D (tipo SIFT)

```cpp
vector<vector<double>> datos_128d;
mt19937 gen(789);
uniform_real_distribution<double> dist(0.0, 255.0);

for(int i = 0; i < 15; i++) {
    vector<double> punto(128);
    for(int j = 0; j < 128; j++) {
        punto[j] = dist(gen);
    }
    datos_128d.push_back(punto);
}

DBfsh indice_128d(128, 1, 2.5, 1, 999);
indice_128d.insertar(datos_128d);

vector<double> query_128d = datos_128d[8];
vector<double> vecino = indice_128d.C_ANN(query_128d, 2.5);
```

### Salida esperada

```
DB-LSH inicializado:
  Dimensión original: 50D
  Dimensión proyectada: 10D (R*-tree 10D)
  Tablas hash: 1
  C = 2, w0 = 16, t = 1
  Semilla: 456

Insertando 12 puntos de 50D...
Proyecciones generadas (primeros 5):
  Punto[0] 50D -> Hash 10D: [0.0107, 0.8406, 0.7435, ...]
  Punto[1] 50D -> Hash 10D: [0.4698, -0.3523, 0.2097, ...]
  ...

c-ANN Query (Algorithm 2)
Query 50D: [1.449, 1.059, 0.037, ...]

--- (r,c)-NN Query ---
Hash G(q) = [-1.199, 1.113, -1.271, 0.005, 0.631, ...]

Tabla 1:
  Window W(G(q), w=16):
    [-9.199, 6.801] × ... × [-8.632, 7.368]
  Puntos encontrados en ventana: 12
  Punto 1: id=0           👈 ¡Usa ID del R*-tree directamente!
    ||q, o|| = 12.536, cr = 2
  ...

RESULTADO:
  Vecino encontrado 50D (primeras 5 dims): [-1.329, 1.131, 0.631, ...]
  Distancia euclidiana: 13.825
```

## 🧪 Funcionamiento detallado de los algoritmos

### **Algorithm 1: (r,c)-NN Query** (del paper)

**Objetivo**: Encontrar un punto `o` tal que `dist(q, o) ≤ c·r`

**Pasos de implementación:**

1. **Proyectar query**: Calcular `G(q) = (h₁(q), h₂(q))` donde:
   ```cpp
   h₁(q) = a₁ · q = Σ(a₁[j] * q[j])
   h₂(q) = a₂ · q = Σ(a₂[j] * q[j])
   ```

2. **Definir ventana**: `W(G(q), w₀·r) = [h₁(q) - w/2, h₁(q) + w/2] × [h₂(q) - w/2, h₂(q) + w/2]`
   - Ancho: `w = w₀ · r` donde `w₀ = 4c²`
   - Centrada en el hash del query

3. **Buscar candidatos**: Ejecutar window query en R*-tree:
   ```cpp
   resultados = indice.windowQuery(x_min, y_min, x_max, y_max);
   ```

4. **Verificar distancias**: Para cada candidato encontrado:
   ```cpp
   int id = res.second;              // ID del R*-tree
   punto_original = datos[id];       // Recuperar original
   dist = distanciaEuclidiana(query, punto_original);
   ```

5. **Condiciones de parada** (del paper):
   - **Éxito**: Si `dist(q, o) ≤ c·r` → retornar `o`
   - **Límite**: Si `cnt = 2tL + 1` → retornar `o` (evita búsqueda infinita)
   - **Fallo**: Si no quedan candidatos → retornar `∅`

### **Algorithm 2: c-ANN Query** (del paper)

**Objetivo**: Encontrar vecino c-aproximado sin conocer la distancia de antemano

**Estrategia de expansión dinámica:**

```
Iteración 1: r = 1    → ventana pequeña
Iteración 2: r = 1.5  → ventana 1.5× más grande (si C=1.5)
Iteración 3: r = 2.25 → ventana 2.25× más grande
...
```

**Flujo completo:**

```
Query q = (6, 6) con c = 1.5

Iter 1: r=1.0, w=9.0  → W(G(q), 9)   → No encuentra → expandir
Iter 2: r=1.5, w=13.5 → W(G(q), 13.5) → Encuentra punto a dist=1.41
        ✓ 1.41 ≤ 1.5·1.5 = 2.25 → RETORNAR punto
```

**Por qué funciona** (garantía del paper):
- Si existe un punto a distancia `d`, eventualmente `r ≥ d/c`
- Entonces la ventana será suficientemente grande para encontrarlo
- La expansión geométrica (`r ← c·r`) garantiza convergencia logarítmica

## 📊 Parámetros del sistema

| Parámetro | Descripción | Valor típico |
|-----------|-------------|--------------||
| **D** | Dimensión original | 20, 50, 128, 700, ... |
| **K** | Dimensión proyectada (fijo) | **10** (para R*-tree 10D) |
| **L** | Tablas hash | 1 (⏳ pendiente L>1) |
| **C** | Factor aproximación | 1.5 - 3.0 |
| **t** | Parámetro tolerancia | 1 |
| **w₀** | Ancho ventana base | 4·C² = 9.0 (C=1.5) |
| **seed** | Semilla aleatoria | Cualquier unsigned int |

## 🔧 Detalles de implementación

### Clase DBfsh (V3)

```cpp
class DBfsh {
private:
    int D;                           // Dimensión original (20, 50, 128, 700, ...)
    int K;                           // Dimensión hash (siempre 10)
    int L;                           // Tablas hash (actualmente 1)
    double C, w0;                    // Parámetros LSH
    int t;                           // Tolerancia
    unsigned seed;                   // Semilla reproducible
    
    vector<vector<double>> a;        // Matriz K×D de proyección (10×D)
    vector<vector<double>> datos;    // Datos originales (índice = id)
    RStarTreeIndex indice;           // R*-tree 10D
    
    void generarFuncionesHash();     // Genera 10 vectores N(0,1)
    array<double,10> funcionHash(const vector<double>& punto);
    
public:
    DBfsh(int dim, int L_, double C_, int t_, unsigned seed_ = 42);
    void insertar(const vector<vector<double>>& datos);
    vector<double> RC_NN(const vector<double>& query, double r, double c);
    vector<double> C_ANN(const vector<double>& query, double c);
};
```

### Flujo de proyección N → 10D

```
Punto original (N dimensiones: 20D, 50D, 128D, 700D, ...)
         ↓
    h₁ = a₁ · p = Σ(a₁[j] * p[j])  → escalar
    h₂ = a₂ · p = Σ(a₂[j] * p[j])  → escalar
    ...
    h₁₀ = a₁₀ · p = Σ(a₁₀[j] * p[j]) → escalar
         ↓
   Hash 10D: [h₁, h₂, ..., h₁₀]
         ↓
    Insertar en R*-tree 10D con ID
```

### **Relación con el paper DB-LSH**

Este proyecto sigue fielmente el paper con estas correspondencias:

| Concepto del paper | Implementación en código |
|-------------------|--------------------------||
| `G_i(q)` | `funcionHash(query)` → `[h₁, h₂, ..., h₁₀]` |
| `W(G_i(q), w₀·r)` | `windowQuery(mins, maxs)` con arrays 10D |
| Buckets L | R*-tree único (⏳ pendiente L>1) |
| Condición `cnt = 2tL+1` | Contador de candidatos inspeccionados |
| Verificación `||q,o||` | `distanciaEuclidiana(query, punto_original)` |
| Expansión `r ← cr` | Loop en `C_ANN` multiplicando r |
| Ecuación (8) ventana | `[h_k(q) - w/2, h_k(q) + w/2]` para k=1..10 |

### **Extensiones más allá del paper**

✅ **N-dimensional**: El paper usa dimensión fija, aquí es configurable  
✅ **K=10 proyecciones**: Mejor discriminación que K=2  
✅ **Vectores aleatorios**: Generados con `N(0,1)` normalizado (reproducible)  
✅ **Optimización O(1)**: Mapeo ID→datos (no mencionado en paper, mejora práctica)  
✅ **Arrays para ventanas**: Interfaz limpia `windowQuery(mins, maxs)`  
✅ **Comentarios del paper**: Código anotado con Algorithm 1/2 líneas  
⏳ **Multi-L pendiente**: L>1 tablas hash independientes (siguiente versión)

```cpp
// Optimización de implementación (no del paper)
void insertar(const vector<vector<double>>& datos_input) {
    datos = datos_input;  // Guardar originales separados
    
    for (size_t i = 0; i < datos.size(); i++) {
        auto hash = funcionHash(datos[i]);
        int id = static_cast<int>(i);        // ID = índice
        indice.insertPrueba(id, hash);       // R*-tree guarda (hash, id)
    }
}

// En RC_NN (del paper, pero optimizado):
int id = res.second;                         // R*-tree da ID (O(1))
const vector<double>& punto = datos[id];     // Acceso directo (O(1))
// En lugar de búsqueda lineal en hash_to_original (O(N))
```

## 🎯 Mejoras futuras

### ⏳ **Pendiente prioritario (V4)**

- ⬜ **Multi-L (L>1)**: Múltiples tablas hash independientes según paper
  - Cada tabla con sus propias K funciones hash G_i
  - Iteración sobre L índices R*-tree separados
  - Mejora garantías probabilísticas del paper

- ⬜ **Carga desde CSV**: Lector robusto para datasets reales
  - Soporte para archivos grandes (60,000+ puntos, 700+ dimensiones)
  - Parsing eficiente con `fstream`
  - Validación de datos y manejo de errores
  - Ejemplo: `datos = cargarCSV("dataset_700d.csv");`

###  (✅ = implementado en V3)

- ✅ **a aleatorios**: N(0,1) normalizados con semilla fija
- ✅ **Separar datos**: Vector `datos` separado del índice (O(1) acceso)
- ✅ **N-dimensional**: Soporta cualquier dimensión de entrada
- ✅ **K=10 dimensiones**: Mejor discriminación que K=2
- ✅ **Arrays ventanas**: `windowQuery(mins, maxs)` limpio
- ⬜ **Métricas**: `recall@k`, `overall ratio`, #candidatos, tiempo


## 📈 Casos de uso

### Datasets típicos

| Dataset | Dimensiones | Tamaño | Uso |
|---------|-------------|--------|-----|
| **Iris** | 4D | 150 | Clasificación |
| **MNIST** | 784D | 60,000 | Dígitos |
| **SIFT** | 128D | 1M - 1B | Imágenes |
| **GloVe** | 50D - 300D | 400K | Word embeddings |
| **Deep1B** | 96D | 1B | Deep features |

### Configuración recomendada

```cpp
// Para SIFT (128D → 10D)
DBfsh sift_index(128, 1, 2.5, 1, 42);

// Para GloVe-200 (200D → 10D)
DBfsh glove_index(200, 1, 2.0, 1, 999);

// Para datasets 700D (alta dimensionalidad)
DBfsh high_dim_index(700, 1, 3.0, 1, 123);

// Para datos densos medianos
DBfsh medium_index(50, 1, 1.5, 1, 456);
```

## 📖 Referencias y conceptos

### **Paper original DB-LSH**

Este proyecto implementa los conceptos de:
- **Dynamic window queries** en lugar de buckets hash estáticos
- **Expansión adaptativa** del radio de búsqueda (`r ← c·r`)
- **Condición de parada dual**: distancia exacta ó límite de candidatos
- **Integración con índices espaciales** (R*-tree en lugar de hash tables)

### **Conceptos clave**

- **LSH (Locality-Sensitive Hashing)**: Proyecta puntos similares a hashes similares
  - Propiedad: Si `dist(p, q)` pequeña → `Pr[h(p) = h(q)]` alta
  - Funciones: `h_i(p) = ⌊(a_i · p + b) / w⌋` (en el paper clásico)
  
- **DB-LSH**: Variante database-friendly
  - NO usa buckets fijos → usa **window queries dinámicas**
  - NO requiere hash tables → usa **índices espaciales existentes**
  - Ajusta ventana según necesidad (expansión `r`)

- **c-ANN**: c-Approximate Nearest Neighbor
  - Encuentra punto `o` tal que `dist(q, o) ≤ c · dist(q, NN)`
  - Donde `NN` es el vecino más cercano real
  - Trade-off: precisión vs velocidad (c=1.5 → acepta 50% error)

- **R*-tree**: Índice espacial multidimensional de Boost.Geometry
  - Organiza puntos en rectángulos MBR (Minimum Bounding Rectangle)
  - Window query: `O(log n + k)` donde k = resultados
  - Permite búsquedas rectangulares eficientes

### **Garantías teóricas del paper**

1. **Probabilidad de éxito**: Con alta probabilidad encuentra vecino c-aproximado
2. **Tiempo**: O(n^ρ log n) donde ρ < 1 depende de c
3. **Espacio**: O(n) para almacenar proyecciones + índice
4. **Aproximación**: Garantiza factor c (configurable)

## 🤝 Contribuciones

Proyecto educativo para Estructuras de Datos Avanzadas. Sugerencias bienvenidas.

## 📄 Licencia

Proyecto académico - Universidad XYZ

---

**Autor**: Manuel J. Silva  
**Versión**: 3.0 (K=10 dimensiones, preparado para L>1 y CSV)  
**Fecha**: Noviembre 2025  
**Tecnologías**: C++17, Boost.Geometry R*-tree 10D, LSH, STL

**Estado actual:**
- ✅ Implementación fiel al paper DB-LSH (Algorithm 1 y 2)
- ✅ Soporte N-dimensional → 10D con K=10 funciones hash
- ✅ Window queries 10D con arrays
- ⏳ Pendiente: Multi-L (L>1) y carga CSV para datasets reales (700D+)
