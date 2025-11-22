#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include "R_star.h"
#include <array> // Agregar include para std::array

using namespace std;
using namespace std::chrono;

// Cargar MNIST desde CSV
vector<vector<double>> cargarMNIST(const string& archivo, int dim_esperada, int limite = -1) {
    vector<vector<double>> dataset;
    ifstream file(archivo);
    string linea;
    
    if(!file.is_open()) {
        cerr << "Error: No se pudo abrir " << archivo << endl;
        return dataset;
    }
    
    cout << "Cargando " << archivo << "..." << endl;
    getline(file, linea); // Skip header
    
    while(getline(file, linea)) {
        if(limite > 0 && dataset.size() >= (size_t)limite) break;
        if(linea.empty()) continue;
        
        vector<double> punto;
        stringstream ss(linea);
        string val;
        
        getline(ss, val, ','); // Skip label
        
        while(getline(ss, val, ',')) {
            try {
                punto.push_back(stod(val));
            } catch(...) {
                continue;
            }
        }
        
        if(punto.size() == (size_t)dim_esperada) {
            dataset.push_back(punto);
        }
    }
    
    file.close();
    cout << "Cargados " << dataset.size() << " puntos de dim=" << dim_esperada << endl;
    return dataset;
}

// Calcular distancia euclidiana
double distanciaEuclidiana(const vector<double>& p1, const vector<double>& p2) {
    double sum = 0.0;
    for(size_t i = 0; i < p1.size(); i++) {
        double diff = p1[i] - p2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

// Validar calidad del resultado
double validarCalidad(const vector<vector<double>>& dataset, 
                      const vector<double>& query, 
                      double dist_encontrada) {
    double min_dist = numeric_limits<double>::max();
    
    for(const auto& punto : dataset) {
        double d = distanciaEuclidiana(query, punto);
        min_dist = min(min_dist, d);
    }
    
    double ratio = dist_encontrada / max(min_dist, 1e-9);
    
    cout << "   Dist DB-LSH: " << dist_encontrada 
         << " | Exacta: " << min_dist 
         << " | Ratio: " << ratio << "x" << endl;
    
    return ratio;
}

// Clase DB-LSH según el paper
class DBLSH {
private:
    int D;              // Dimensión original
    int K;              // Dimensión proyectada
    int L;              // Número de tablas hash
    int t;              // Parámetro de balance
    double C;           // Approximation ratio c
    double w0;          // Ancho de bucket inicial
    double p1, p2;      // Probabilidades de colisión
    unsigned seed;
    
    vector<RStarTreeIndex> indices;
    vector<vector<double>> datos;
    vector<vector<vector<double>>> familias_hash; // [L][K][D]
    
    // Función de distribución normal estándar
    double phi(double x) {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    }
    
    // Probabilidad de colisión según Ecuación (4) del paper
    double calcularP(double tau, double w) {
        // p(τ; w) = ∫_{-w/(2τ)}^{w/(2τ)} f(t) dt
        double limite = w / (2.0 * tau);
        return 2.0 * phi(limite) - 1.0;
    }
    
    void generarFuncionesHash() {
        familias_hash.resize(L);
        mt19937 gen(seed);
        normal_distribution<double> dist(0.0, 1.0);
        
        for(int l = 0; l < L; l++) {
            familias_hash[l].resize(K);
            for(int k = 0; k < K; k++) {
                familias_hash[l][k].resize(D);
                for(int d = 0; d < D; d++) {
                    familias_hash[l][k][d] = dist(gen);
                }
            }
        }
    }
    
    // Función hash según Ecuación (3)
    vector<double> funcionHash(const vector<double>& punto, int l_idx) {
        vector<double> resultado(K);
        for(int k = 0; k < K; k++) {
            resultado[k] = 0.0;
            for(int d = 0; d < D; d++) {
                resultado[k] += familias_hash[l_idx][k][d] * punto[d];
            }
        }
        return resultado;
    }

public:
    DBLSH(int dim, int K_fijo, int L_fijo, double c, int t_param, unsigned seed_val = 42) 
        : D(dim), K(K_fijo), L(L_fijo), t(t_param), C(c), w0(4*c*c), seed(seed_val) {
        
        // Calcular probabilidades según Sección III-C
        p1 = calcularP(1.0, w0);
        p2 = calcularP(C, w0);
        
        cout << "\n=== Parámetros DB-LSH ===" << endl;
        cout << "D (dim original): " << D << endl;
        cout << "K (dim proyectada): " << K << " [FIJO]" << endl;
        cout << "L (num tablas): " << L << " [FIJO]" << endl;
        cout << "c (approx ratio): " << C << endl;
        cout << "t (balance): " << t << endl;
        cout << "w0 = 4c²: " << w0 << endl;
        cout << "p1 = p(1; w0): " << p1 << endl;
        cout << "p2 = p(c; w0): " << p2 << endl;
    }
    
    void insertar(const vector<vector<double>>& datos_entrada) {
        datos = datos_entrada;
        int n = datos.size();
        
        // Calcular ρ* para información (no se usa para ajustar K y L)
        double rho_star = log(1.0 / p1) / log(1.0 / p2);
        
        cout << "\n=== Información del Dataset ===" << endl;
        cout << "n (datos): " << n << endl;
        cout << "ρ* teórico: " << rho_star << endl;
        cout << "K (FIJO): " << K << endl;
        cout << "L (FIJO): " << L << endl;
        cout << "Límite candidatos: " << (2*t*L + 1) << endl;
        
        // Inicializar R*-Trees
        indices.clear();
        indices.resize(L);
        
        // Generar funciones hash
        generarFuncionesHash();
        
        // Bulk Loading (crítico para eficiencia)
        cout << "\nRealizando Bulk Loading..." << endl;
        
        for(int l = 0; l < L; l++) {
            vector<pair<vector<double>, int>> entradas;
            entradas.reserve(n);
            
            for(int i = 0; i < n; i++) {
                vector<double> hash_val = funcionHash(datos[i], l);
                entradas.push_back({hash_val, i});
            }
            
            // Sort-Tile-Recursive (STR) packing implícito
            sort(entradas.begin(), entradas.end(), 
                 [](const auto& a, const auto& b) {
                     return a.first[0] < b.first[0];
                 });
            
            indices[l].bulkInsert(entradas);
        }
        
        cout << "Indexación completada." << endl;
    }
    
    // Algoritmo 1: (r,c)-NN Query
    pair<vector<double>, bool> RC_NN(const vector<double>& query, double r) {
        int cnt = 0;
        double w_r = w0 * r;
        double threshold = w_r / 2.0;
        
        vector<double> ultimo_candidato;
        int limite_candidatos = 2 * t * L + 1;
        
        for(int l = 0; l < L; l++) {
            vector<double> hash_query = funcionHash(query, l);
            
            // Construir ventana W(G_i(q), w0*r) según Ecuación (8)
            vector<double> mins(K), maxs(K);
            for(int k = 0; k < K; k++) {
                mins[k] = hash_query[k] - threshold;
                maxs[k] = hash_query[k] + threshold;
            }
            
            // Window Query en R*-Tree
            vector<RStarTreeIndex::Value> resultados = indices[l].windowQuery(
                *reinterpret_cast<const array<double, 10>*>(mins.data()),
                *reinterpret_cast<const array<double, 10>*>(maxs.data())
            );
            
            for(const auto& res : resultados) {
                int id = res.second; // Corregir acceso al ID
                ultimo_candidato = datos[id];
                double dist = distanciaEuclidiana(query, datos[id]);
                cnt++;
                
                // Condición de terminación 1 (Línea 6 del Algoritmo 1)
                if(dist <= C * r) {
                    return {datos[id], true};
                }
                
                // Condición de terminación 2 (Línea 6 del Algoritmo 1)
                if(cnt >= limite_candidatos) {
                    return {ultimo_candidato, true};
                }
            }
        }
        
        // Caso 2 de Definición 2: no hay punto dentro de c*r
        return {vector<double>(), false};
    }
    
    // Algoritmo 2: c-ANN Query
    pair<vector<double>, double> C_ANN(const vector<double>& query) {
        double r = 1.0;
        
        // Calcular límite superior de búsqueda
        double max_dist = 0.0;
        for(size_t i = 0; i < min(size_t(100), datos.size()); i++) {
            double d = distanciaEuclidiana(query, datos[i]);
            max_dist = max(max_dist, d);
        }
        double r_max = max_dist * 10.0;
        
        // Búsqueda incremental con r = 1, c, c², ...
        while(r < r_max) {
            auto [resultado, encontrado] = RC_NN(query, r);
            
            if(encontrado && !resultado.empty()) {
                double dist = distanciaEuclidiana(query, resultado);
                return {resultado, dist};
            }
            
            r *= C;  // Incrementar radio exponencialmente
        }
        
        return {vector<double>(), -1.0};
    }
    
    void printStats() {
        cout << "\n=== Estadísticas DB-LSH ===" << endl;
        cout << "Puntos indexados: " << datos.size() << endl;
        cout << "K (dim proyectada): " << K << endl;
        cout << "L (num R*-Trees): " << L << endl;
        if(!indices.empty()) {
            indices[0].printStats();
        }
    }
};

// MAIN
int main() {
    // Configuración FIJA
    const string CSV_FILE = "fashion-mnist_test.csv";
    const int DIM = 784;
    const int LIMITE_DATOS = 5000;
    const int NUM_QUERIES = 20;
    const double C = 1.5;
    const unsigned SEED = 42;
    
    // PARÁMETROS FIJOS (según paper para datasets medianos)
    const int K_FIJO = 10;  // Dimensión proyectada
    const int L_FIJO = 5;   // Número de tablas hash
    
    // 1. Cargar datos
    cout << "=== CARGANDO DATOS ===" << endl;
    vector<vector<double>> full_dataset = cargarMNIST(CSV_FILE, DIM, LIMITE_DATOS);
    if(full_dataset.empty()) {
        cerr << "Error: No se pudieron cargar datos" << endl;
        return 1;
    }
    
    // Shuffle
    shuffle(full_dataset.begin(), full_dataset.end(), default_random_engine(SEED));
    
    // Separar en dataset de indexación y queries
    vector<vector<double>> dataset_index, dataset_queries;
    for(size_t i = 0; i < full_dataset.size() - NUM_QUERIES; i++)
        dataset_index.push_back(full_dataset[i]);
    for(size_t i = full_dataset.size() - NUM_QUERIES; i < full_dataset.size(); i++)
        dataset_queries.push_back(full_dataset[i]);
    
    // 2. Probar diferentes valores de T (balance eficiencia/calidad)
    ofstream resultados("resultados_dblsh.csv");
    resultados << "T,TiempoPromedio_ms,RatioPromedio,Recall" << endl;
    
    vector<int> valores_T = {5, 10, 20, 50};
    
    for(int T : valores_T) {
        cout << "\n\n========================================" << endl;
        cout << ">>> PROBANDO T = " << T << " <<<" << endl;
        cout << "========================================" << endl;
        
        // Crear índice con K y L FIJOS
        DBLSH db(DIM, K_FIJO, L_FIJO, C, T, SEED);
        db.insertar(dataset_index);
        db.printStats();
        
        // Ejecutar queries
        double tiempo_total_ms = 0.0;
        double ratio_total = 0.0;
        int exitos = 0;
        
        cout << "\nEjecutando " << NUM_QUERIES << " queries..." << endl;
        for(size_t q = 0; q < dataset_queries.size(); q++) {
            cout << "\n[Query " << (q+1) << "/" << NUM_QUERIES << "]" << endl;
            
            auto start = high_resolution_clock::now();
            auto [resultado, dist] = db.C_ANN(dataset_queries[q]);
            auto stop = high_resolution_clock::now();
            
            double tiempo_ms = duration_cast<microseconds>(stop - start).count() / 1000.0;
            tiempo_total_ms += tiempo_ms;
            
            if(!resultado.empty() && dist >= 0) {
                double ratio = validarCalidad(dataset_index, dataset_queries[q], dist);
                ratio_total += ratio;
                exitos++;
            } else {
                cout << "   No encontrado" << endl;
            }
        }
        
        // Resultados
        double tiempo_prom = tiempo_total_ms / NUM_QUERIES;
        double ratio_prom = (exitos > 0) ? (ratio_total / exitos) : 0.0;
        double recall = (double)exitos / NUM_QUERIES;
        
        cout << "\n--- RESUMEN T=" << T << " ---" << endl;
        cout << "Tiempo promedio: " << tiempo_prom << " ms" << endl;
        cout << "Ratio promedio: " << ratio_prom << "x" << endl;
        cout << "Recall: " << recall << " (" << exitos << "/" << NUM_QUERIES << ")" << endl;
        
        resultados << T << "," << tiempo_prom << "," << ratio_prom << "," << recall << endl;
    }
    
    resultados.close();
    
    cout << "\n========================================" << endl;
    cout << "Resultados guardados en 'resultados_dblsh.csv'" << endl;
    cout << "========================================" << endl;
    
    return 0;
}
