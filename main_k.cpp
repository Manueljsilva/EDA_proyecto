#include <iostream>
#include <fstream>
#include <sstream>
#include "R_star2.h"
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <random>
#include <set>
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
        double R_min; // Radio mínimo inicial
        double t;
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
        // Ground Truth: Encontrar k vecinos más cercanos reales (fuerza bruta)
        vector<pair<int, double>> encontrarKVecinosReales(const vector<double>& query, int k) {
            // Calcular distancias a todos los puntos
            vector<pair<double, int>> distancias;
            distancias.reserve(datos.size());
            
            for(size_t i = 0; i < datos.size(); i++) {
                double dist = distanciaEuclidiana(query, datos[i]);
                distancias.push_back({dist, i});
            }
            
            // Ordenar por distancia (parcial sort hasta k)
            partial_sort(distancias.begin(), 
                        distancias.begin() + min(k, (int)distancias.size()),
                        distancias.end());
            
            // Retornar los k primeros como {id, distancia}
            vector<pair<int, double>> resultado;
            int num_vecinos = min(k, (int)distancias.size());
            for(int i = 0; i < num_vecinos; i++) {
                resultado.push_back({distancias[i].second, distancias[i].first});
            }
            
            return resultado;
        }
        
        // Constructor con parámetros del paper DB-LSH (según código original)
        // D: dimensión original, L: número de tablas hash, C: approximation ratio
        // R_min: radio inicial mínimo, beta: proporción máxima de accesos (0-1)
        // w0 = R_min * 4C² según código original del paper
        DBLSH(int dim, int L_, double C_, double R_min_, double t_, unsigned seed_ = 42) 
            : D(dim), L(L_), C(C_), R_min(R_min_), t(t_), seed(seed_) {
            w0 = R_min * 4.0 * C * C;  // Fórmula del código original
            indices.resize(L);
            generarFuncionesHash();
            
            cout << "DB-LSH inicializado (según implementación original):" << endl;
            cout << "  Dimensión original: " << D << "D" << endl;
            cout << "  Dimensión proyectada: " << K << "D (R*-tree " << K << "D)" << endl;
            cout << "  Tablas hash: " << L << endl;
            cout << "  C = " << C << ", R_min = " << R_min << ", t = " << t << endl;
            cout << "  w0 = " << w0 << " (R_min * 4C²)" << endl;
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
                }
                // Bulk-loading: construir R*-tree de una sola vez (más eficiente según paper DB-LSH)
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
        
        
        // Algorithm 1 (modificado): (r,c)-NN Query para k vecinos
        // Input: q (query point), r (query radius), c (approximation ratio), k (num neighbors), T (límite de accesos)
        // Output: Lista de hasta k puntos con {id, punto, distancia}
        vector<tuple<int, vector<double>, double>> RC_NN_K(const vector<double>& query, double r, double c, int k, int T){
            vector<tuple<int, vector<double>, double>> candidatos; // {id, punto, distancia}
            set<int> ids_visitados; // Evitar duplicados entre tablas
            int cnt = 0;
            
            // Para cada tabla i = 1 to L
            for(int i = 0; i < L; i++){
                array<double, K> hash_query = funcionHash(query, i);
                double w_r = w0 * r;
                double threshold = w_r / 2.0;
                
                array<double, K> mins, maxs;
                for(size_t j = 0; j < K; j++) {
                    mins[j] = hash_query[j] - threshold;
                    maxs[j] = hash_query[j] + threshold;
                }
                
                vector<typename RStarTreeIndex<K>::Value> resultados = indices[i].windowQuery(mins, maxs);

                for(const auto& res : resultados) {
                    int id = res.second;
                    
                    // Evitar duplicados
                    if(ids_visitados.count(id)) continue;
                    ids_visitados.insert(id);
                    
                    const vector<double>& punto_original = datos[id];
                    double dist = distanciaEuclidiana(query, punto_original);
                    cnt++;
                    
                    // Agregar si dist ≤ cr
                    if(dist <= c * r) {
                        candidatos.push_back({id, punto_original, dist});
                        if((int)candidatos.size() >= k) {
                            return candidatos; // Terminación temprana si ya tenemos k
                        }
                    }
                    
                    // Terminación por límite de accesos T (según código original)
                    if(cnt >= T) {
                        return candidatos;
                    }
                }
            }
            return candidatos;
        }
        
        // Algorithm 2 (modificado según código original): c-ANN Query para k vecinos
        // Input: q (query point), c (approximation ratio), k (num neighbors)
        // Output: Lista de k puntos con {id, punto, distancia}
        vector<tuple<int, vector<double>, double>> C_ANN_K(const vector<double>& query, double c, int k){
            // Calcular parámetro t adaptativo según tamaño N (código original)
            // int N = datos.size();
            //double t = 1.0;
            // if (N < 70000) {
            //     t = 200.0;
            // } else if (N >= 70000 && N < 500000) {
            //     t = 1000.0;
            // } else if (N >= 500000 && N < 2000000) {
            //     t = 2000.0;
            // } else if (N >= 2000000 && N < 2000000000) {
            //     t = 20000.0;
            // } else {
            //     t = 20000.0;
            // }
            // t *= 2.0;
            
            int T = 2*t*L + k;
            double r = R_min;

            
            // Acumular candidatos entre iteraciones con IDs
            vector<tuple<int, vector<double>, double>> acumulados;
            set<int> ids_usados; // Para evitar duplicados usando IDs reales
            
            // int rounds = 0;
            // const int MAX_ROUNDS = 30;  // Límite de 30 rondas (código original)
            
            while(true){
                // rounds++;
                auto nuevos = RC_NN_K(query, r, c, k, T);  // r = init_w/w0
                
                // Agregar nuevos candidatos evitando duplicados
                for(const auto& candidato : nuevos) {
                    int id = get<0>(candidato);
                    
                    if(ids_usados.find(id) == ids_usados.end()) {
                        ids_usados.insert(id);
                        acumulados.push_back(candidato);
                    }
                }
                
                if((int)acumulados.size() >= k) {
                    // Ordenar por distancia y retornar los k mejores
                    sort(acumulados.begin(), acumulados.end(), 
                         [](const auto& a, const auto& b) { return get<2>(a) < get<2>(b); });
                    
                    acumulados.resize(min(k, (int)acumulados.size()));
                    return acumulados;
                }
                
                // Expandir ventana w para siguiente iteración (código original)
                r *= c;
            }
            
            // Si no se encontraron k vecinos después de MAX_ROUNDS, retornar lo acumulado
            sort(acumulados.begin(), acumulados.end(),
                 [](const auto& a, const auto& b) { return get<2>(a) < get<2>(b); });
            return acumulados;
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

// ============ MÉTRICAS DEL PAPER ============

// Distancia euclidiana (helper fuera de clase)
double distEuclidiana(const vector<double>& p1, const vector<double>& p2) {
    double sum = 0.0;
    for(size_t i = 0; i < p1.size(); i++) {
        double diff = p1[i] - p2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

// Overall Ratio para k-NN (Ecuación 11 del paper)
// Promedio de ratios individuales: ratio_i = dist(q, o_i) / dist(q, o*_i)
double calcularOverallRatioKNN(const vector<tuple<int, vector<double>, double>>& vecinos_dblsh,
                                const vector<pair<int, double>>& vecinos_reales) {
    if(vecinos_dblsh.empty() || vecinos_reales.empty()) return 1.0;
    
    double sum_ratios = 0.0;
    int n = min(vecinos_dblsh.size(), vecinos_reales.size());
    
    for(int i = 0; i < n; i++) {
        double dist_dblsh = get<2>(vecinos_dblsh[i]);
        double dist_real = vecinos_reales[i].second;
        if(dist_real > 1e-9) {
            sum_ratios += dist_dblsh / dist_real;
        }
    }
    
    return sum_ratios / n;
}

// Recall para k-NN (Ecuación 12 del paper) - VERSIÓN CORRECTA
// Recall = |R ∩ R*| / k
// donde R = IDs devueltos por DB-LSH, R* = IDs de vecinos reales
double calcularRecallKNN(const vector<tuple<int, vector<double>, double>>& vecinos_dblsh,
                         const vector<pair<int, double>>& vecinos_reales) {
    if(vecinos_dblsh.empty() || vecinos_reales.empty()) return 0.0;
    
    // Extraer IDs de DB-LSH (R)
    set<int> ids_dblsh;
    for(const auto& v : vecinos_dblsh) {
        ids_dblsh.insert(get<0>(v));
    }
    
    // Extraer IDs de vecinos reales (R*)
    set<int> ids_reales;
    for(const auto& v : vecinos_reales) {
        ids_reales.insert(v.first);
    }
    
    // Calcular intersección |R ∩ R*|
    int interseccion = 0;
    for(int id : ids_dblsh) {
        if(ids_reales.count(id)) {
            interseccion++;
        }
    }
    
    int k = ids_reales.size();
    return (double)interseccion / k;
}

#define MAIN_C 1.01
#define MAIN_t 500
#define MAIN_K 68
#define MAIN_L 18

int main(){
    std::filesystem::create_directories("results");

    cout << "============================================================" << endl;
    cout << "DB-LSH: Evaluación con k=50 Queries" << endl;
    cout << "============================================================" << endl;

    cout << "Cargando Fashion-MNIST...\n";
    vector<vector<double>> full_dataset = loadDataset("fashion_mnist.csv", 60000);
    cout << "Filas cargadas: " << full_dataset.size() << endl;

    const int D = 784;
    const double C = MAIN_C;
    const int t = MAIN_t;
    const int K = MAIN_K;
    const int L = MAIN_L;
    // Separar queries del dataset de indexación
    const int K_QUERIES = 50;
    const double R_MIN = 1;  // Radio inicial mínimo (parámetro típico del paper)
    
    // Shuffle dataset
    mt19937 gen(42);
    shuffle(full_dataset.begin(), full_dataset.end(), gen);
    
    // Tomar primeras K_QUERIES como queries
    vector<vector<double>> queries;
    for(int i = 0; i < K_QUERIES; i++) {
        queries.push_back(full_dataset[i]);
    }
    
    // Resto para indexar
    vector<vector<double>> dataset_index;
    for(size_t i = K_QUERIES; i < full_dataset.size(); i++) {
        dataset_index.push_back(full_dataset[i]);
    }
    
    cout << "Queries: " << queries.size() << endl;
    cout << "Dataset indexado: " << dataset_index.size() << endl;

    // Construir índice DB-LSH (con parámetros del código original)
    DBLSH<K> indice(D, L, C, R_MIN, t, 42); // D=784, L=5 , C=1.5, R_min=0.3, beta=0.1, seed=42 
    indice.insertar(dataset_index);
    indice.imprimir();

    // ============ EJECUTAR k-NN QUERIES ============
    // Probar con diferentes valores de k (como en la Figura 5 del paper)
    vector<int> k_values = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    
    cout << "\n" << string(60, '=') << endl;
    cout << "EJECUTANDO k-NN QUERIES (1 query, diferentes k)" << endl;
    cout << string(60, '=') << endl;
    
    // Usar primera query como ejemplo
    const auto& query = queries[0];
    
    // Abrir archivo CSV para guardar resultados
    ofstream csv_file("results/knn_results.csv");
    csv_file << "k,recall,overall_ratio\n";
    
    cout << "\nResultados para diferentes valores de k:\n" << endl;
    cout << "k\tRecall\t\tOverall Ratio" << endl;
    cout << string(60, '-') << endl;
    
    for(int k : k_values) {
        cout << "\rProcesando k=" << k << "..." << flush;
        
        // 1. Buscar k vecinos con DB-LSH (ahora retorna tuplas con IDs)
        auto vecinos_dblsh = indice.C_ANN_K(query, C, k);
        
        // 2. Ground truth: k vecinos reales
        auto vecinos_reales = indice.encontrarKVecinosReales(query, k);
        
        // 3. Calcular métricas
        double recall = calcularRecallKNN(vecinos_dblsh, vecinos_reales);
        double ratio = calcularOverallRatioKNN(vecinos_dblsh, vecinos_reales);
        
        // Guardar en CSV
        csv_file << k << "," << recall << "," << ratio << "\n";
    }
    
    cout << "\r" << string(60, ' ') << "\r";  // Limpiar línea de progreso
    csv_file.close();
    
    // Releer y mostrar resultados
    ifstream csv_read("results/knn_results.csv");
    string line;
    getline(csv_read, line);  // Skip header
    
    while(getline(csv_read, line)) {
        stringstream ss(line);
        string k_str, recall_str, ratio_str;
        getline(ss, k_str, ',');
        getline(ss, recall_str, ',');
        getline(ss, ratio_str, ',');
        cout << k_str << "\t" << recall_str << "\t\t" << ratio_str << endl;
    }
    csv_read.close();
    
    cout << "\n[Resultados guardados en results/knn_results.csv]" << endl;
    
    cout << "\n" << string(60, '=') << endl;
    cout << "\n[Interpretación]" << endl;
    cout << "- Recall: |R ∩ R*| / k (fracción de IDs que coinciden)" << endl;
    cout << "- Overall Ratio: promedio de dist_dblsh/dist_real para k vecinos" << endl;
    cout << "- Valores cercanos a 1.0 (ratio) y 1.0 (recall) indican mejor calidad" << endl;
    cout << string(60, '=') << endl;

    return 0;
}
