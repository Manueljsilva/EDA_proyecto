# DB-LSH: Implementación Fiel al Paper Original

**Locality-Sensitive Hashing with Query-based Dynamic Bucketing (ICDE 2022)**

## 📋 ¿Qué es DB-LSH?

**DB-LSH** es una variante de LSH diseñada para integrarse con estructuras de indexación espacial existentes (como R*-tree) en lugar de usar tablas hash tradicionales.

---

## 📖 Descripción

Este proyecto implementa **fielmente el algoritmo DB-LSH** según el pseudocódigo publicado en el paper ICDE 2022.

**DB-LSH** combina Locality-Sensitive Hashing con índices espaciales R*-tree para búsquedas de vecinos más cercanos aproximados (c-ANN) en alta dimensión.

---

## 🎯 Características

- ✅ **Implementación fiel al paper**: Algoritmos 1 y 2 exactos
- ✅ **Fashion-MNIST**: 60,000 imágenes 784D
- ✅ **Proyección LSH**: 784D → K-D (configurable)
- ✅ **R*-tree espacial**: Índices multidimensionales con bulk-loading
- ✅ **Métricas del paper**: Recall (Eq. 12), Overall Ratio (Eq. 11)
- ✅ **Análisis comparativo**: Paper vs GitHub original
- ✅ **Documentación completa**: Diferencias y mejoras identificadas

---

## 🚀 Inicio Rápido

### Requisitos

```bash
# Ubuntu/Debian
sudo apt-get install g++ make libboost-all-dev python3 python3-pip

# Python (visualización)
pip3 install pandas matplotlib numpy
```

### Compilación y Ejecución

```bash
# Compilar
make clean && make all

# Experimento k-NN (k=1,10,20,...,100)
make run-k

# Experimento varying n (Fig. 5, 6, 7 del paper)
make run-grafico
```

---

## 📊 Experimentos Disponibles

### 1. k-NN Benchmark (`main_k.cpp`)

Evalúa calidad de resultados para diferentes valores de k:

```bash
make run-k
```

**Configuración:**
- k ∈ {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}
- **Parámetros usados:**
  - `C = 1.01` (approximation ratio)
  - `t = 500` (parámetro de límite de accesos)
  - `K = 68` (dimensión proyectada para R*-tree)
  - `L = 18` (número de tablas hash)
  - `n = 60,000` (tamaño del dataset)

**Salida:**
- `results/knn_results.csv`: Recall y Ratio por k

---

### 2. Varying n (`main_grafico.cpp`)

Replica experimentos del paper variando tamaño del dataset:

```bash
make run-grafico
```

**Configuración:**
- n ∈ {0.2, 0.4, 0.6, 0.8, 1.0} (proporción del dataset)
- k = 50 vecinos
- 50 queries promediadas
- **Parámetros usados:**
  - `C = 1.5` (approximation ratio)
  - `t = 8000` (parámetro de límite de accesos)
  - `K = 83` (dimensión proyectada para R*-tree)
  - `L = 2` (número de tablas hash)

**Salida:**
- `results/varying_n_results.csv`: n, recall, ratio, tiempo promedio

---




### Resultados Experimentales

#### Experimento k-NN (Configuración: C=1.01, t=500, K=68, L=18)

| k | Recall (%) | Overall Ratio |
|---|------------|---------------|
| 10 | ~3-5% | ~1.01x |
| 50 | ~2-4% | ~1.01x |
| 100 | ~1-3% | ~1.01x |

*Nota: Ejecutar `make run-k` para obtener resultados detallados*

#### Experimento Varying n (Configuración: C=1.5, t=8000, K=83, L=2, k=50)

| n (proporción) | Dataset Size | Recall (%) | Ratio | Time (ms) |
|----------------|--------------|------------|-------|----------|
| 0.2 | 12,000 | ~35% | ~1.14x | ~0.9 |
| 0.6 | 36,000 | ~30% | ~1.16x | ~1.1 |
| 1.0 | 60,000 | ~25% | ~1.18x | ~1.2 |

*Nota: Ejecutar `make run-grafico` para obtener resultados detallados*

**Observaciones:**
- Configuraciones diferentes producen trade-offs distintos entre recall y eficiencia
- C más bajo (1.01) mejora ratio pero reduce recall
- C más alto (1.5) mejora recall pero aumenta ratio
- El recall disminuye con datasets más grandes (comportamiento esperado en LSH)



---

## 📂 Estructura del Proyecto

```
EDA_proyecto/
├── R_star2.h                    # Implementación R*-tree con Boost.Geometry
├── main_k.cpp                   # Experimento k-NN benchmark
├── main_grafico.cpp             # Experimento varying n
├── main.cpp                     # Testing sintético
├── Makefile                     # Compilación y ejecución
├── fashion_mnist.csv            # Dataset Fashion-MNIST (60k imágenes)
├── calculate_parameters.py      # Script para calcular parámetros óptimos
├── pgm_to_png.py               # Utilidad de conversión de imágenes
├── README.md                    # Este archivo
├── bin/                         # Ejecutables compilados
│   ├── main
│   ├── main_k
│   └── main_grafico
├── obj/                         # Archivos objeto (.o)
└── results/                     # Resultados experimentales
    ├── knn_results.csv         # Resultados k-NN benchmark
    └── varying_n_results.csv   # Resultados varying n
```

---

## ⚙️ Configuración de Parámetros

### Modificar Parámetros (main_k.cpp, main_grafico.cpp)

Los parámetros se definen con `#define` al final de cada archivo:

```cpp
// main_k.cpp - Ejemplo configuración óptima encontrada
#define MAIN_C 1.01      // Approximation ratio
#define MAIN_t 500       // Parámetro t del paper
#define MAIN_K 68        // Dimensión proyectada (R*-tree)
#define MAIN_L 18        // Número de tablas hash

// main_grafico.cpp - Ejemplo configuración
#define MAIN_C 1.5
#define MAIN_t 8000
#define MAIN_K 83
#define MAIN_L 2
```

### Parámetros Fijos en Constructor

```cpp
const int D = 784;           // Dimensión Fashion-MNIST
const double R_MIN = 1.0;    // Radio mínimo inicial
const unsigned seed = 42;    // Reproducibilidad
```

### Fórmulas del Paper

```cpp
// Límite de accesos (Algorithm 1)
T = 2 * t * L + 1

// Ancho de ventana inicial
w₀ = R_min * 4 * C²

// Expansión de radio (Algorithm 2)
r_new = r_old * C
```

---

## 🎯 Algoritmos Implementados

### Algorithm 1: (r,c)-NN Query

```cpp
// Input: query q, radio r, approx ratio c, k neighbors, límite T
vector<tuple<int, vector<double>, double>> RC_NN_K(q, r, c, k, T) {
    candidatos = [];
    cnt = 0;
    
    for (i = 0; i < L; i++) {  // Para cada tabla hash
        W = windowQuery(hash(q), w₀·r);  // Búsqueda en R*-tree
        
        for (o in W) {
            cnt++;
            dist = ||q - o||;
            
            // PAPER: Filtro por distancia (línea 10)
            if (dist <= c·r) {
                candidatos.push_back({id, punto, dist});
            }
            
            if (cnt >= T || |candidatos| >= k) {
                return candidatos;
            }
        }
    }
    return candidatos;
}
```

### Algorithm 2: c-ANN Query

```cpp
// Input: query q, approx ratio c, k neighbors
vector<tuple<int, vector<double>, double>> C_ANN_K(q, c, k) {
    r = R_min;
    T = 2*t*L + 1;  // PAPER: Límite fijo
    acumulados = [];
    
    while (true) {  // PAPER: Sin límite de rondas
        nuevos = RC_NN_K(q, r, c, k, T);
        acumular(nuevos);  // Evitar duplicados
        
        if (|acumulados| >= k) {
            return top_k(acumulados);  // Ordenar por distancia
        }
        
        r = r * c;  // Expandir radio
    }
}
```

**Nota:** Esta implementación sigue el pseudocódigo exacto. Ver comentarios en código para diferencias con GitHub.

---

## 📈 Métricas Implementadas

### 1. Recall (Ecuación 12)

$$\text{Recall} = \frac{|R \cap R^*|}{k}$$

- **R**: Conjunto de IDs devueltos por DB-LSH
- **R***: Conjunto de IDs de k vecinos reales (ground truth)
- **Implementación**: Intersección de sets

```cpp
set<int> ids_dblsh = {IDs de resultado DB-LSH};
set<int> ids_reales = {IDs de ground truth};
recall = intersection(ids_dblsh, ids_reales).size() / k;
```

### 2. Overall Ratio (Ecuación 11)

$$\text{Overall Ratio} = \frac{1}{k} \sum_{i=1}^{k} \frac{\text{dist}(q, o_i)}{\text{dist}(q, o_i^*)}$$

- Promedio de ratios individuales de distancias
- **Ideal**: 1.0 (exacto)
- **Aceptable**: ≤ c (approximation ratio)

```cpp
ratio = mean(dist_dblsh[i] / dist_real[i] for i in [0..k-1]);
```

### 3. Query Time

```cpp
auto start = high_resolution_clock::now();
vecinos = indice.C_ANN_K(query, C, k);
auto end = high_resolution_clock::now();
time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
```

---

## 🔧 Detalles de Implementación

### Clase Principal

```cpp
template <size_t K>  // K = dimensión proyectada (R*-tree)
class DBLSH {
private:
    int D;                              // Dimensión original (784)
    int L;                              // Número de tablas hash
    double C, w0, R_min, t;
    vector<RStarTreeIndex<K>> indices;  // L R*-trees de K dimensiones
    vector<vector<double>> datos;       // Datos originales D-dim
    vector<vector<vector<double>>> a;   // Proyecciones L×K×D
    
public:
    DBLSH(int dim, int L, double C, double R_min, double t, unsigned seed);
    void insertar(const vector<vector<double>>& datos);
    
    // Algorithm 1
    vector<tuple<int, vector<double>, double>> 
        RC_NN_K(query, r, c, k, T);
    
    // Algorithm 2  
    vector<tuple<int, vector<double>, double>> 
        C_ANN_K(query, c, k);
    
    // Ground truth (fuerza bruta)
    vector<pair<int, double>> 
        encontrarKVecinosReales(query, k);
};
```

### Proyección LSH

```cpp
// Hash functions gaussianas N(0,1)
for (i = 0; i < L; i++)
    for (j = 0; j < K; j++)
        for (d = 0; d < D; d++)
            a[i][j][d] ~ N(0,1)

// Proyectar punto D-dim → K-dim
hash_result[j] = Σ(a[tabla][j][d] * punto[d]) for d in [0..D-1]
```

### R*-tree Bulk-Loading

```cpp
// Construir índice eficientemente
for (tabla = 0; tabla < L; tabla++) {
    proyecciones = [];
    for (punto in datos) {
        hash = funcionHash(punto, tabla);  // D → K dimensiones
        proyecciones.push_back({hash, id});
    }
    indices[tabla].bulkLoad(proyecciones);  // Construcción bottom-up
}
```



---

## 🔧 Troubleshooting

### Problemas Comunes

**Error: `No se pudo abrir fashion_mnist.csv`**
```bash
# Verificar que el archivo existe en el directorio raíz
ls -lh fashion_mnist.csv
```

**Compilación falla con error de Boost**
```bash
# Instalar Boost en Ubuntu/Debian
sudo apt-get install libboost-all-dev

# Verificar versión (requiere 1.70+)
dpkg -l | grep libboost
```

**Resultados con recall 0% o muy bajo**
- Verificar parámetros: C, t, K, L deben ser apropiados para el dataset
- Probar con `calculate_parameters.py` para encontrar configuración óptima
- El recall bajo puede indicar que `t` es muy pequeño o `C` muy estricto

**Tiempo de ejecución muy largo**
- Reducir el tamaño del dataset en `loadDataset()`
- Ajustar parámetros: reducir `L` o `t`
- Compilar con optimizaciones: `make clean && make all`

**Memoria insuficiente**
```bash
# Reducir dataset en main_k.cpp o main_grafico.cpp
# Cambiar: loadDataset("fashion_mnist.csv", 60000)
# Por:     loadDataset("fashion_mnist.csv", 30000)
```

---

## 📚 Documentación Complementaria

### Archivos Importantes

1. **Código fuente principal**
   - `R_star2.h`: Implementación del índice R*-tree
   - `main_k.cpp`: Experimento k-NN con parámetros optimizados
   - `main_grafico.cpp`: Experimento varying n del paper
   - `main.cpp`: Testing con dataset sintético

2. **Comentarios en código**
   - Headers explicativos en main_k.cpp y main_grafico.cpp
   - Documentación inline de algoritmos
   - Referencias a ecuaciones del paper (Eq. 11, 12)

3. **Resultados experimentales**
   - `results/knn_results.csv`: k, recall, ratio
   - `results/varying_n_results.csv`: n, recall, ratio, time

4. **Scripts de utilidad**
   - `calculate_parameters.py`: Optimización de parámetros C, t, K, L
   - `pgm_to_png.py`: Conversión de imágenes PGM a PNG

---

## 🛠️ Comandos Make

```bash
# Compilar todos
make all

# Ejecutables individuales
make main_k          # k-NN benchmark
make main_grafico    # Varying n
make main2           # Testing sintético

# Ejecutar
make run-k           # k-NN experiments
make run-grafico     # Varying n experiments
make run-test        # main2 (dataset sintético)

# Utilidades
make clean           # Limpiar binarios
make rebuild         # Limpiar y recompilar
```

---

## 📖 Referencias

### Paper Original

- **Título:** DB-LSH: Locality-Sensitive Hashing with Query-based Dynamic Bucketing
- **Autores:** Yao Tian, Xi Zhao, Xiaofang Zhou
- **Conferencia:** IEEE ICDE 2022
- **DOI:** 10.1109/ICDE53745.2022.00264
- **GitHub:** https://github.com/Jacyhust/DB-LSH

### Dataset

- **Fashion-MNIST** by Zalando Research
- **Tamaño:** 60,000 training images + 10,000 test
- **Formato:** 28×28 grayscale (784 features)
- **Clases:** 10 (T-shirt, Trouser, Pullover, Dress, Coat, Sandal, Shirt, Sneaker, Bag, Ankle boot)

### Tecnologías

- **C++17:** std::filesystem, structured bindings, templates
- **Boost 1.70+:** Boost.Geometry R*-tree
- **Python 3:** pandas, matplotlib, numpy


---

## 👨‍💻 Autores

**Silva Anampa, Manuel Jesus**  
**Campoverde San Martín, Yacira Nicol**  
**Bracamonte Toguchi, Mikel Dan**  
**Garcia Calle, Renato**  
Proyecto de Estructuras de Datos Avanzadas  
Curso: EDA  
Fecha: Noviembre 2024

---
