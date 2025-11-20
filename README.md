# DB-LSH N-Dimensional (V2) con R*-Tree

Implementación de **DB-LSH** (Database-friendly Locality-Sensitive Hashing) siguiendo el paper original, con extensiones para **N dimensiones**. Utiliza **window queries dinámicas** sobre un R*-tree en lugar de buckets estáticos, permitiendo búsquedas aproximadas de vecinos más cercanos (c-ANN) eficientes.

## 📋 ¿Qué es DB-LSH?

**DB-LSH** es una variante de LSH diseñada para integrarse con estructuras de indexación espacial existentes (como R*-tree) en lugar de usar tablas hash tradicionales.

### Diferencias clave con LSH clásico:

| Característica | LSH clásico | DB-LSH (este proyecto) |
|----------------|-------------|------------------------|
| **Almacenamiento** | Buckets hash estáticos | R*-tree dinámico |
| **Búsqueda** | Lookup directo en bucket | Window query expansiva |
| **Colisiones** | Todos los puntos en bucket | Ventana `W(G(q), w₀·r)` |
| **Expansión** | Probar múltiples tablas L | Expandir radio `r ← c·r` |
| **Indexación** | Hash tables | Índice espacial (R*-tree) |

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

## 🔬 Implementación en este proyecto

Esta versión extiende DB-LSH con soporte **N-dimensional** y optimizaciones de rendimiento:

- **Proyecciones LSH**: Funciones hash `h_i(p) = a_i · p` con vectores `a_i ~ N(0,1)` normalizados
- **N → 2D**: Reduce cualquier dimensión a 2D para indexar en R*-tree 2D de Boost
- **Window queries**: Usa `windowQuery(x_min, y_min, x_max, y_max)` del R*-tree
- **Verificación final**: Calcula distancia euclidiana real en espacio original (N-D)
- **Acceso O(1)**: Optimización usando ID del R*-tree (no en el paper, mejora práctica)

## ✨ Mejoras de V2 vs V1

### ✅ **1. Soporte N-dimensional dinámico**

```cpp
// V1: Solo 2D fijo
DBfsh indice_2d(2, 1, 2, 1.5, 1);

// V2: Cualquier dimensión
DBfsh indice_5d(5, 1, 1.5, 1);    // 5D → 2D
DBfsh indice_10d(10, 1, 1.5, 1);  // 10D → 2D
DBfsh indice_128d(128, 1, 1.5, 1); // 128D → 2D (SIFT)
```

### ✅ **2. Vectores hash aleatorios N(0,1)**

**Antes (V1):** Hardcodeado
```cpp
vector<vector<double>> a = {{0.6, 0.8}, {0.3, -0.9}};  // Fijo
```

**Ahora (V2):** Generación aleatoria normalizada
```cpp
// Genera K vectores de D dimensiones ~ N(0,1), normalizados
void generarFuncionesHash() {
    mt19937 gen(seed);  // Reproducible
    normal_distribution<double> dist(0.0, 1.0);
    // ... normalización ||a[i]|| = 1
}
```

### ✅ **3. Optimización O(1) para recuperación de datos**

**❌ Antes (V1):** O(N) búsqueda lineal
```cpp
// Mapeo ineficiente
vector<pair<tuple<double,double>, vector<double>>> hash_to_original;

// Búsqueda O(N) comparando hashes con tolerancia
for(const auto& [hash, original] : hash_to_original) {
    if(abs(get<0>(hash) - get<0>(hash_candidato)) < 0.0001) {
        punto_original = original;  // ¡60,000 comparaciones!
    }
}
```

**✅ Ahora (V2):** O(1) acceso directo
```cpp
// Vector simple (índice = id)
vector<vector<double>> datos;

// Inserción: usa índice como ID
int id = static_cast<int>(i);
indice.insertPrueba(id, hash_punto);

// Recuperación: O(1) usando ID del R*-tree
int id = res.second;  // R*-tree da el ID directamente
const vector<double>& punto_original = datos[id];  // ¡Acceso directo!
```

**Comparación de rendimiento:**

| Operación | V1 (hash_to_original) | V2 (vector + ID) |
|-----------|----------------------|------------------|
| **Recuperación** | O(N) búsqueda lineal | **O(1) acceso directo** |
| **60,000 puntos** | 60,000 comparaciones | 1 acceso indexado |
| **Comparación doubles** | Insegura (tolerancia 0.0001) | No necesaria |
| **Memoria** | 2× (duplica hash+original) | 1× (solo originales) |
| **Escalabilidad** | ❌ Empeora con N | ✅ Constante siempre |

### ✅ **4. Estructura de datos genérica**

```cpp
// V1: tuple<double, double> (solo 2D)
vector<tuple<double, double>> datos_2d;

// V2: vector<double> (N dimensiones)
vector<vector<double>> datos_nd = {
    {1.0, 2.0, 3.0, 4.0, 5.0},     // 5D
    {1.0, 2.0, ..., 128.0}         // 128D
};
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
V_2/
├── main.cpp           # Implementación N-dimensional optimizada
├── R_star.h           # Interfaz del R*-tree
├── R_star.cpp         # Implementación del R*-tree
├── Makefile           # Sistema de compilación
├── README.md          # Este archivo
├── bin/               # Ejecutables (generado)
└── obj/               # Archivos objeto (generado)
```

## 🔍 Ejemplos de uso

### Ejemplo 1: 2D → 2D (caso base)

```cpp
vector<vector<double>> datos_2d = {
    {1.0, 1.0}, {2.0, 2.0}, {4.0, 2.0},
    {5.0, 5.0}, {7.0, 8.0}
};

DBfsh indice_2d(2, 1, 1.5, 1, 42);  // dim=2, L=1, C=1.5, t=1, seed=42
indice_2d.insertar(datos_2d);

vector<double> query = {6.0, 6.0};
vector<double> vecino = indice_2d.C_ANN(query, 1.5);
```

### Ejemplo 2: 5D → 2D

```cpp
vector<vector<double>> datos_5d = {
    {1.0, 2.0, 3.0, 4.0, 5.0},
    {2.0, 3.0, 4.0, 5.0, 6.0},
    {5.0, 5.0, 5.0, 5.0, 5.0}
};

DBfsh indice_5d(5, 1, 1.5, 1, 123);  // 5 dimensiones
indice_5d.insertar(datos_5d);

vector<double> query_5d = {1.2, 2.1, 3.2, 4.1, 5.1};
vector<double> vecino = indice_5d.C_ANN(query_5d, 1.5);
```

### Ejemplo 3: 10D → 2D

```cpp
// Generar puntos 10D
vector<vector<double>> datos_10d;
for(int i = 0; i < 100; i++) {
    vector<double> punto(10);
    for(int j = 0; j < 10; j++) {
        punto[j] = rand() / double(RAND_MAX) * 100.0;
    }
    datos_10d.push_back(punto);
}

DBfsh indice_10d(10, 1, 1.5, 1, 999);
indice_10d.insertar(datos_10d);

vector<double> query_10d(10, 5.0);  // Query de 10 dimensiones
vector<double> vecino = indice_10d.C_ANN(query_10d, 2.0);
```

### Salida esperada

```
DB-LSH inicializado:
  Dimensión original: 5D
  Dimensión proyectada: 2D (fijo para R*-tree)
  Tablas hash: 1
  C = 1.5, w0 = 9, t = 1
  Semilla: 123

Insertando 5 puntos de 5D...
Proyecciones generadas (primeros 5):
  Punto[0] 5D -> Hash: (6.36023, -3.26102)
  Punto[1] 5D -> Hash: (8.06503, -3.85684)
  ...

c-ANN Query
Query 5D: [1.2, 2.1, 3.2, 4.1, 5.1]

--- (r,c)-NN Query ---
Ventana W(G(q), w_0·r = 9): [2.05, 11.05] x [-7.89, 1.11]
  Puntos encontrados en ventana: 4
  Punto 1: id=0           👈 ¡Usa ID del R*-tree directamente!
    dist(q, o) = 0.331662
  ✓ Condición ||q,o|| ≤ cr cumplida

RESULTADO:
  Vecino encontrado 5D: [1, 2, 3, 4, 5]
  Distancia euclidiana: 0.331662
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
|-----------|-------------|--------------|
| **D** | Dimensión original | 2, 5, 10, 128, ... |
| **K** | Dimensión proyectada (fijo) | 2 (para R*-tree 2D) |
| **L** | Tablas hash | 1 (simplificado) |
| **C** | Factor aproximación | 1.5 - 2.0 |
| **t** | Parámetro tolerancia | 1 |
| **w₀** | Ancho ventana base | 4·C² = 9.0 (C=1.5) |
| **seed** | Semilla aleatoria | Cualquier unsigned int |

## 🔧 Detalles de implementación

### Clase DBfsh (V2)

```cpp
class DBfsh {
private:
    int D;                           // Dimensión original (configurable)
    int K;                           // Dimensión hash (siempre 2)
    int L;                           // Tablas hash
    double C, w0;                    // Parámetros LSH
    int t;                           // Tolerancia
    unsigned seed;                   // Semilla reproducible
    
    vector<vector<double>> a;        // Matriz K×D de proyección
    vector<vector<double>> datos;    // Datos originales (índice = id)
    RStarTreeIndex indice;           // R*-tree
    
    void generarFuncionesHash();     // Genera vectores N(0,1)
    tuple<double,double> funcionHash(const vector<double>& punto);
    
public:
    DBfsh(int dim, int L_, double C_, int t_, unsigned seed_ = 42);
    void insertar(const vector<vector<double>>& datos);
    vector<double> RC_NN(const vector<double>& query, double r, double c);
    vector<double> C_ANN(const vector<double>& query, double c);
};
```

### Flujo de proyección N → 2D

```
Punto original (N dimensiones)
         ↓
    h₁ = a₁ · p = Σ(a₁[j] * p[j])  → escalar
    h₂ = a₂ · p = Σ(a₂[j] * p[j])  → escalar
         ↓
   Hash 2D: (h₁, h₂)
         ↓
    Insertar en R*-tree con ID
```

### **Relación con el paper DB-LSH**

Este proyecto sigue fielmente el paper con estas correspondencias:

| Concepto del paper | Implementación en código |
|-------------------|--------------------------|
| `G_i(q)` | `funcionHash(query)` → `(h₁, h₂)` |
| `W(G_i(q), w₀·r)` | `windowQuery(x_min, y_min, x_max, y_max)` |
| Buckets L | R*-tree único (simplificado L=1) |
| Condición `cnt = 2tL+1` | Contador de candidatos inspeccionados |
| Verificación `\|\|q,o\|\|` | `distanciaEuclidiana(query, punto_original)` |
| Expansión `r ← cr` | Loop en `C_ANN` multiplicando r |

### **Extensiones más allá del paper**

✅ **N-dimensional**: El paper usa dimensión fija, aquí es configurable  
✅ **Vectores aleatorios**: Generados con `N(0,1)` normalizado (reproducible)  
✅ **Optimización O(1)**: Mapeo ID→datos (no mencionado en paper, mejora práctica)

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

### Corto plazo (✅ = implementado en V2)

- ✅ **a aleatorios**: N(0,1) normalizados con semilla fija
- ✅ **Separar datos**: Vector `datos` separado del índice (O(1) acceso)
- ✅ **N-dimensional**: Soporta cualquier dimensión de entrada
- ⬜ **Deduplicación de candidatos**: `unordered_set<int>` entre iteraciones
- ⬜ **CSV (fstream)**: Lector robusto para datasets reales
- ⬜ **Métricas**: `recall@k`, `overall ratio`, #candidatos, tiempo

### Medio plazo

- ⬜ **k-ANN**: Con `priority_queue` (max-heap tamaño k)
- ⬜ **Multi-L**: Múltiples tablas hash (L > 1)
- ⬜ **Parámetros por CLI**: `--w0 --c --r0 --L --seed`
- ⬜ **Tests unitarios**: Proyección, ventana, verificación

### Largo plazo

- ⬜ **Multi-probe**: Ventanas adyacentes priorizadas
- ⬜ **Persistencia**: Guardar/cargar índice serializado
- ⬜ **Batch queries**: Procesar múltiples queries eficientemente
- ⬜ **Benchmarks**: SIFT-1M, GloVe-200, Deep1B

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
// Para SIFT (128D → 2D)
DBfsh sift_index(128, 1, 1.5, 1, 42);

// Para GloVe-200 (200D → 2D)
DBfsh glove_index(200, 1, 2.0, 1, 999);

// Para datos densos pequeños
DBfsh small_index(10, 1, 1.2, 1, 123);
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

**Autor**: Manuel J. Simpson  
**Versión**: 2.0 (Optimizado N-dimensional)  
**Fecha**: Noviembre 2025  
**Tecnologías**: C++17, Boost.Geometry, R*-tree, LSH, STL
