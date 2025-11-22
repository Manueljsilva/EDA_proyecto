#include <iostream>
#include <fstream>
#include <sstream>
#include "R_star2.h"
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <random>
#include <filesystem>

using namespace std;
// clase  
template <size_t K> // Número de funciones hash = dimensión proyectada (ahora K=10)
class DBLSH {
    private:
        int D;      // Dimensión original (ej: 2, 10, 128, 700)
        int L;      // Número de tablas hash
        double C;   // Constante de aproximación
        double w0;  // Ancho de ventana base
        int t;      // Parámetro t (entero positivo)
        unsigned seed; // Semilla para reproducibilidad
        
        vector<RStarTreeIndex<K>> indices;
        
        // Almacena los puntos originales (índice del vector = id del punto)
        vector<vector<double>> datos;

        // Matriz de proyección: K filas × D columnas
        // a[i][j] = coeficiente de la función hash i para la dimensión j
        vector<vector<vector<double>>> a;

        // Generar funciones hash aleatorias
        void generarFuncionesHash() {
            a.resize(L);
            mt19937 gen(seed);
            normal_distribution<double> dist(0.0, 1.0);
            for(int i = 0; i < L; i++) {
                a[i].resize(K);
                
                for(size_t j = 0; j < K; j++) {
                    a[i][j].resize(D);
                    
                    // Generar vector aleatorio N(0,1)
                    for(int k = 0; k < D; k++) {
                        a[i][j][k] = dist(gen);
                    }
                }
            }
        }
        
        // Proyectar punto N-dimensional → K-dimensional (K=10)
        array<double, K> funcionHash(const vector<double>& punto, int tabla) {
            if(static_cast<int>(punto.size()) != D) {
                throw runtime_error("Punto debe tener " + to_string(D) + " dimensiones");
            }
            
            array<double, K> hash_result;
            
            // h_i(p) = a[i] · p (producto punto) para cada función hash
            for(size_t i = 0; i < K; i++) {
                hash_result[i] = 0.0;
                for(int j = 0; j < D; j++) {
                    hash_result[i] += a[tabla][i][j] * punto[j];
                }
            }
            
            return hash_result;
        }
        
        // Distancia euclidiana en espacio original N-dimensional
        double distanciaEuclidiana(const vector<double>& p1, const vector<double>& p2) {
            double sum = 0.0;
            for(size_t i = 0; i < p1.size(); i++) {
                double diff = p1[i] - p2[i];
                sum += diff * diff;
            }
            return sqrt(sum);
        }


    public:
        // Constructor con parámetros del paper DB-LSH
        // D: dimensión original, L: número de tablas hash, C: approximation ratio
        // t: parámetro para termination condition (cnt = 2tL+1)
        // w0 = 4C² según sugerencia del paper para (1,c,p1,p2)-sensitive hash family
        DBLSH(int dim, int L_, double C_, int t_, unsigned seed_ = 42) 
            : D(dim), L(L_), C(C_), w0(4*C_*C_), t(t_), seed(seed_) {
            indices.resize(L);
            generarFuncionesHash();
            
            cout << "DB-LSH inicializado:" << endl;
            cout << "  Dimensión original: " << D << "D" << endl;
            cout << "  Dimensión proyectada: " << K << "D (R*-tree " << K << "D)" << endl;
            cout << "  Tablas hash: " << L << endl;
            cout << "  C = " << C << ", w0 = " << w0 << ", t = " << t << endl;
            cout << "  Semilla: " << seed << endl;
        }

        void insertar(const vector<vector<double>>& datos_input){
            // Guardar datos originales
            datos = datos_input;
            
            cout << "\nIndexando " << datos.size() << " puntos de " << D << "D..." << endl;
            cout << "Usando bulk-loading (paper DB-LSH)" << endl;
            
            // Proyectar TODOS los puntos primero (preparación para bulk-loading)
            vector<vector<pair<array<double, K>, int>>> proyecciones;
            proyecciones.resize(L);
            // Proyectar y guardar en R*-tree con ID = índice del vector
            for(int i = 0; i < L; i++) {
                proyecciones[i].reserve(datos.size());
                for (size_t j = 0; j < datos.size(); j++){
                    array<double, K> hash_punto = funcionHash(datos[j], i);
                    int id = static_cast<int>(j);  // ID = índice en vector datos
                    proyecciones[i].push_back({hash_punto, id});
                    indices[i].insertPrueba(id, hash_punto);
                }
                // Bulk-loading: construir R*-tree de una sola vez (más eficiente)
                indices[i].bulkLoad(proyecciones[i]);
            }
                        
            cout << "Proyecciones generadas (primeros 5):" << endl;
            for (size_t i = 0; i < min((size_t)5, datos.size()); i++) {
                auto hash = funcionHash(datos[i], 0);
                cout << "  Punto[" << i << "] " << D << "D -> Hash "<< K << "D: [";
                for(size_t k = 0; k < K; k++) {
                    cout << hash[k];
                    if(k < K-1) cout << ", ";
                }
                cout << "]" << endl;
            }
            if(datos.size() > 5) {
                cout << "  ... (" << (datos.size() - 5) << " más)" << endl;
            }
            cout << endl;
        }
        void imprimir(){
            for(int i = 0; i < L; i++) {
                indices[i].printStats();
            }
        }
        
        
        // Algorithm 1: (r,c)-NN Query
        // Input: q (query point), r (query radius), c (approximation ratio), t (positive integer)
        // Output: A point o or ∅
        vector<double> RC_NN(const vector<double>& query, double r, double c){
            int cnt = 0;
            
            cout << "\n--- (r,c)-NN Query ---" << endl;
            cout << "Query " << D << "D: [";
            for(size_t i = 0; i < min((size_t)5, query.size()); i++) {
                cout << query[i];
                if(i < min((size_t)4, query.size()-1)) cout << ", ";
            }
            if(query.size() > 5) cout << " ...";
            cout << "]" << endl;
            cout << "r = " << r << ", c = " << c << ", t = " << t << ", L = " << L << endl;
            
            // Para cada tabla i = 1 to L (Algorithm 1, line 2)
            for(int i = 0; i < L; i++){
                cout << "\nTabla " << (i+1) << ":" << endl;
                
                // Compute G_i(q) (Algorithm 1, line 3)
                // G(q) = (h_1(q), h_2(q), ..., h_K(q)) proyección K-dimensional
                array<double, K> hash_query = funcionHash(query, i);
                cout << "Hash G(q) = [";
                for(size_t k = 0; k < K; k++) {
                    cout << hash_query[k];
                    if(k < K-1) cout << ", ";
                }
                cout << "]" << endl;
                
                // Window query W(G_i(q), w_0·r) (Algorithm 1, line 4)
                // Ecuación (8): W(G(q), w) = [h_1(q) - w/2, h_1(q) + w/2] × ... × [h_K(q) - w/2, h_K(q) + w/2]
                double w_r = w0 * r;  // w = w_0 · r
                double threshold = w_r / 2.0;  // w/2 para cada dimensión
                
                // Construir hiper-rectángulo K-dimensional
                array<double, K> mins, maxs;
                for(size_t k = 0; k < K; k++) {
                    mins[k] = hash_query[k] - threshold;  // h_k(q) - w/2
                    maxs[k] = hash_query[k] + threshold;  // h_k(q) + w/2
                }
                
                cout << "  Window W(G(q), w=" << w_r << "):" << endl;
                cout << "    [" << mins[0] << ", " << maxs[0] << "] × ... × [" << mins[K - 1] << ", " << maxs[K - 1] << "]" << endl;
                
                // Window Query en el R*-tree 10D (ahora con arrays)
                vector<typename RStarTreeIndex<K>::Value> resultados = indices[i].windowQuery(mins, maxs);
                cout << "  Puntos encontrados en ventana: " << resultados.size() << endl;
                
                // Procesar cada punto encontrado (Algorithm 1, line 4-7)
                for(const auto& res : resultados) {
                    // Recuperar ID del punto (el R*-tree retorna pair<Point_hash, id>)
                    int id = res.second;
                    
                    // Recuperar punto ORIGINAL en espacio D-dimensional (O(1))
                    const vector<double>& punto_original = datos[id];
                    
                    cnt = cnt + 1;  // Algorithm 1, line 5
                    cout << "  Punto " << cnt << ": id=" << id << endl;
                    
                    // Verificar distancia en espacio ORIGINAL (Algorithm 1, line 6)
                    // if cnt = 2tL + 1 or ||q, o|| ≤ cr then return o;
                    double dist = distanciaEuclidiana(query, punto_original);
                    cout << "    ||q, o|| = " << dist << ", cr = " << c*r << endl;
                    
                    if(cnt == 2*t*L + 1) {
                        cout << "  ✓ Condición cnt = 2tL+1 = " << (2*t*L + 1) << " alcanzada" << endl;
                        return punto_original;
                    }
                    
                    if(dist <= c * r) {
                        cout << "  ✓ Condición ||q,o|| ≤ cr cumplida" << endl;
                        return punto_original;
                    }
                }
            }
            
            // return ∅ (Algorithm 1, line 8)
            cout << "  ∅ No se encontró punto satisfactorio después de " << cnt << " verificaciones" << endl;
            return vector<double>(D, INFINITY); // Representamos ∅ con infinitos
        }
        
        // Algorithm 2: c-ANN Query (del paper DB-LSH)
        // Input: q (query point), c (approximation ratio)
        // Output: A point o
        vector<double> C_ANN(const vector<double>& query, double c){
            cout << "\n" << string(60, '=') << endl;
            cout << "c-ANN Query (Algorithm 2)" << endl;
            cout << string(60, '=') << endl;
            cout << "Query " << D << "D: [";
            for(size_t i = 0; i < min((size_t)5, query.size()); i++) {
                cout << query[i];
                if(i < min((size_t)4, query.size()-1)) cout << ", ";
            }
            if(query.size() > 5) cout << " ...";
            cout << "]" << endl;
            cout << "c = " << c << endl;
            
            double r = 1.0; // r ← 1 (Algorithm 2, line 1)
            int iteracion = 0;
            
            // while TRUE do (Algorithm 2, line 2)
            while(true){
                iteracion++;
                cout << "\n--- Iteración " << iteracion << " ---" << endl;
                cout << "r = " << r << endl;
                
                // o ← call (r,c)-NN (Algorithm 2, line 3)
                vector<double> o = RC_NN(query, r, c);
                
                // if o ≠ ∅ then return o (Algorithm 2, lines 4-5)
                if(o[0] != INFINITY) {
                    cout << "\n✓ Punto encontrado!" << endl;
                    return o;
                }
                // else r ← cr (Algorithm 2, line 7)
                else {
                    r = c * r;
                    cout << "  No encontrado, expandiendo radio: r ← c·r = " << r << endl;
                }
                
                // Límite de seguridad
                if(r > 1000.0) {
                    cout << "  ⚠ Alcanzado límite de expansión (r > 1000)" << endl;
                    return vector<double>(D, INFINITY);
                }
            }
        }


};

vector<vector<double>> loadDataset(const string& path, size_t max_rows = 5000) {
    ifstream f(path);
    if (!f) throw runtime_error("No se pudo abrir " + path);

    string line;
    getline(f, line); // header
    vector<vector<double>> datos;
    datos.reserve(max_rows);

    while (getline(f, line) && datos.size() < max_rows) {
        stringstream ss(line);
        string cell;

        // salta label
        if (!getline(ss, cell, ',')) continue;

        vector<double> row;
        row.reserve(784);
        while (getline(ss, cell, ',')) {
            row.push_back(stod(cell));
        }
        if (row.size() == 784) datos.push_back(move(row));
    }
    return datos;
}

void guardarPGM(const vector<double>& v, const string& path) {
    if (v.size() != 784) throw runtime_error("Esperaba 784 valores");
    ofstream out(path, ios::binary);
    if (!out) throw runtime_error("No pude abrir " + path);
    out << "P2\n28 28\n255\n"; // formato ASCII PGM
    for (size_t i = 0; i < v.size(); ++i) {
        int val = static_cast<int>(round(clamp(v[i], 0.0, 255.0)));
        out << val << (i % 28 == 27 ? "\n" : " ");
    }
}

int main(){
    std::filesystem::create_directories("results");

    cout << "============================================================" << endl;
    cout << "DB-LSH Dinámico: Reducción de Dimensionalidad N-D → 10D" << endl;
    cout << "============================================================" << endl;

    cout << "Cargando Fashion-MNIST...\n";
    vector<vector<double>> datos = loadDataset("fashion_mnist.csv", 10000);

    cout << "Filas cargadas: " << datos.size() << endl;

    // <K>, D, L, C, T, seed
    DBLSH<10> indice(784, 5, 2.0, 1, 42);
    indice.insertar(datos);
    indice.imprimir();

    // Query: usa la primera imagen como ejemplo con un poco de ruido
    vector<double> query = datos.front();
    for (double& v : query) v += 0.05; // pequeño ruido

    guardarPGM(query, "results/query.pgm");

    auto vecino = indice.C_ANN(query, 2.0);
    // imprime primeras dimensiones
    guardarPGM(vecino, "results/result.pgm");

    return 0;

    /*
    
    // ============= EJEMPLO 1: 20D → 10D =============
    cout << "\n### EJEMPLO 1: Datos 20D → Hash 10D (reducción 20→10) ###\n" << endl;
    
    vector<vector<double>> datos_20d;
    // Generar 8 puntos en 20D con patrones reconocibles
    for(int i = 0; i < 8; i++) {
        vector<double> punto(20);
        for(int j = 0; j < 20; j++) {
            punto[j] = (i + 1) * 0.5 + j * 0.1;  // Valores progresivos
        }
        datos_20d.push_back(punto);
    }

    DBfsh indice_20d(20, 1, 1.5, 1, 42);  // dim=20, L=1, C=1.5, t=1, seed=42
    indice_20d.insertar(datos_20d);
    indice_20d.imprimir();

    // Query similar al punto 3 (índice 2)
    vector<double> query_20d(20);
    for(int j = 0; j < 20; j++) {
        query_20d[j] = 3 * 0.5 + j * 0.1 + 0.05;  // Pequeña perturbación
    }
    
    cout << "\nQuery 20D (primeras 5 dims): [";
    for(int i = 0; i < 5; i++) cout << query_20d[i] << (i < 4 ? ", " : " ...]");
    cout << endl;
    
    vector<double> vecino_20d = indice_20d.C_ANN(query_20d, 1.5);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO:" << endl;
    if(vecino_20d[0] != INFINITY) {
        cout << "  Vecino encontrado 20D (primeras 5 dims): [";
        for(size_t i = 0; i < min((size_t)5, vecino_20d.size()); i++) {
            cout << vecino_20d[i];
            if(i < 4) cout << ", ";
        }
        cout << " ...]" << endl;
        
        // Calcular distancia
        double dist = 0.0;
        for(size_t i = 0; i < query_20d.size(); i++) {
            double diff = query_20d[i] - vecino_20d[i];
            dist += diff * diff;
        }
        dist = sqrt(dist);
        cout << "  Distancia euclidiana: " << dist << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;
    
    
    // ============= EJEMPLO 2: 50D → 10D =============
    cout << "\n\n### EJEMPLO 2: Datos 50D → Hash 10D (reducción 50→10) ###\n" << endl;
    
    vector<vector<double>> datos_50d;
    mt19937 gen(123);
    normal_distribution<double> dist_normal(0.0, 1.0);
    
    // Generar 12 puntos aleatorios en 50D
    for(int i = 0; i < 12; i++) {
        vector<double> punto(50);
        for(int j = 0; j < 50; j++) {
            punto[j] = dist_normal(gen) + i * 0.2;  // Offset por cluster
        }
        datos_50d.push_back(punto);
    }
    
    DBfsh indice_50d(50, 1, 2.0, 1, 456);  // dim=50, L=1, C=2.0, t=1
    indice_50d.insertar(datos_50d);
    indice_50d.imprimir();
    
    // Query similar al punto índice 5
    vector<double> query_50d = datos_50d[5];
    // Agregar ruido
    for(int j = 0; j < 50; j++) {
        query_50d[j] += dist_normal(gen) * 0.1;
    }
    
    cout << "\nBuscando vecino cercano en espacio 50D..." << endl;
    cout << "Query 50D (primeras 5 dims): [";
    for(int i = 0; i < 5; i++) cout << query_50d[i] << (i < 4 ? ", " : " ...]");
    cout << endl;
    
    vector<double> vecino_50d = indice_50d.C_ANN(query_50d, 2.0);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO:" << endl;
    if(vecino_50d[0] != INFINITY) {
        cout << "  Vecino encontrado 50D (primeras 5 dims): [";
        for(size_t i = 0; i < min((size_t)5, vecino_50d.size()); i++) {
            cout << vecino_50d[i];
            if(i < 4) cout << ", ";
        }
        cout << " ...]" << endl;
        
        double dist = 0.0;
        for(size_t i = 0; i < query_50d.size(); i++) {
            double diff = query_50d[i] - vecino_50d[i];
            dist += diff * diff;
        }
        dist = sqrt(dist);
        cout << "  Distancia euclidiana: " << dist << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;
    
    
    // ============= EJEMPLO 3: 128D → 10D (SIFT-like) =============
    cout << "\n\n### EJEMPLO 3: Datos 128D → Hash 10D (reducción 128→10, tipo SIFT) ###\n" << endl;
    
    vector<vector<double>> datos_128d;
    mt19937 gen2(789);
    uniform_real_distribution<double> dist_uniform(0.0, 255.0);
    
    // Generar 15 descriptores tipo SIFT en 128D
    for(int i = 0; i < 15; i++) {
        vector<double> punto(128);
        for(int j = 0; j < 128; j++) {
            punto[j] = dist_uniform(gen2);
        }
        datos_128d.push_back(punto);
    }
    
    DBfsh indice_128d(128, 1, 2.5, 1, 999);  // dim=128, L=1, C=2.5, t=1
    indice_128d.insertar(datos_128d);
    indice_128d.imprimir();
    
    // Query = punto 8 con pequeña variación
    vector<double> query_128d = datos_128d[8];
    for(int j = 0; j < 128; j++) {
        query_128d[j] += (dist_uniform(gen2) - 127.5) * 0.05;  // ±5% ruido
    }
    
    cout << "\nBuscando vecino cercano en espacio 128D (tipo SIFT)..." << endl;
    cout << "Query 128D (primeras 8 dims): [";
    for(int i = 0; i < 8; i++) cout << (int)query_128d[i] << (i < 7 ? ", " : " ...]");
    cout << endl;
    
    vector<double> vecino_128d = indice_128d.C_ANN(query_128d, 2.5);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO:" << endl;
    if(vecino_128d[0] != INFINITY) {
        cout << "  Vecino encontrado 128D (primeras 8 dims): [";
        for(size_t i = 0; i < min((size_t)8, vecino_128d.size()); i++) {
            cout << (int)vecino_128d[i];
            if(i < 7) cout << ", ";
        }
        cout << " ...]" << endl;
        
        double dist = 0.0;
        for(size_t i = 0; i < query_128d.size(); i++) {
            double diff = query_128d[i] - vecino_128d[i];
            dist += diff * diff;
        }
        dist = sqrt(dist);
        cout << "  Distancia euclidiana L2: " << dist << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;
    
    cout << "\n✨ Demostración: Reducción de dimensionalidad 20D/50D/128D → 10D ✨" << endl;
    cout << "💡 Preparado para datasets reales de alta dimensionalidad (SIFT, GloVe, etc.)" << endl;

    return 0;

    */
}
