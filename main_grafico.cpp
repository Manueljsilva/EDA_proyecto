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
#include <chrono>

using namespace std;
using namespace std::chrono;

// Clase DB-LSH (copiada de main_k.cpp)
template <size_t K>
class DBLSH {
    private:
        int D;
        int L;
        double C;
        double w0;
        double R_min;
        double t;
        unsigned seed;
        
        vector<RStarTreeIndex<K>> indices;
        vector<vector<double>> datos;
        vector<vector<vector<double>>> a;

        void generarFuncionesHash() {
            a.resize(L);
            mt19937 gen(seed);
            normal_distribution<double> dist(0.0, 1.0);
            for(int i = 0; i < L; i++) {
                a[i].resize(K);
                for(size_t j = 0; j < K; j++) {
                    a[i][j].resize(D);
                    for(int k = 0; k < D; k++) {
                        a[i][j][k] = dist(gen);
                    }
                }
            }
        }
        
        array<double, K> funcionHash(const vector<double>& punto, int tabla) {
            if(static_cast<int>(punto.size()) != D) {
                throw runtime_error("Punto debe tener " + to_string(D) + " dimensiones");
            }
            
            array<double, K> hash_result;
            for(size_t i = 0; i < K; i++) {
                hash_result[i] = 0.0;
                for(int j = 0; j < D; j++) {
                    hash_result[i] += a[tabla][i][j] * punto[j];
                }
            }
            return hash_result;
        }
        
        double distanciaEuclidiana(const vector<double>& p1, const vector<double>& p2) {
            double sum = 0.0;
            for(size_t i = 0; i < p1.size(); i++) {
                double diff = p1[i] - p2[i];
                sum += diff * diff;
            }
            return sqrt(sum);
        }

    public:
        vector<pair<int, double>> encontrarKVecinosReales(const vector<double>& query, int k) {
            vector<pair<double, int>> distancias;
            distancias.reserve(datos.size());
            
            for(size_t i = 0; i < datos.size(); i++) {
                double dist = distanciaEuclidiana(query, datos[i]);
                distancias.push_back({dist, i});
            }
            
            partial_sort(distancias.begin(), 
                        distancias.begin() + min(k, (int)distancias.size()),
                        distancias.end());
            
            vector<pair<int, double>> resultado;
            int num_vecinos = min(k, (int)distancias.size());
            for(int i = 0; i < num_vecinos; i++) {
                resultado.push_back({distancias[i].second, distancias[i].first});
            }
            return resultado;
        }
        
        DBLSH(int dim, int L_, double C_, double R_min_, double t_, unsigned seed_ = 42) 
            : D(dim), L(L_), C(C_), R_min(R_min_), t(t_), seed(seed_) {
            w0 = R_min * 4.0 * C * C;
            indices.resize(L);
            generarFuncionesHash();
        }

        void insertar(const vector<vector<double>>& datos_input){
            datos = datos_input;
            vector<vector<pair<array<double, K>, int>>> proyecciones;
            proyecciones.resize(L);
            
            for(int i = 0; i < L; i++) {
                proyecciones[i].reserve(datos.size());
                for (size_t j = 0; j < datos.size(); j++){
                    array<double, K> hash_punto = funcionHash(datos[j], i);
                    int id = static_cast<int>(j);
                    proyecciones[i].push_back({hash_punto, id});
                }
                indices[i].bulkLoad(proyecciones[i]);
            }
        }
        
        vector<tuple<int, vector<double>, double>> RC_NN_K(const vector<double>& query, double r, double c, int k, int T){
            vector<tuple<int, vector<double>, double>> candidatos;
            set<int> ids_visitados;
            int cnt = 0;
            
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
                    if(ids_visitados.count(id)) continue;
                    ids_visitados.insert(id);
                    
                    const vector<double>& punto_original = datos[id];
                    double dist = distanciaEuclidiana(query, punto_original);
                    cnt++;
                    
                    if(dist <= c * r) {
                        candidatos.push_back({id, punto_original, dist});
                        if((int)candidatos.size() >= k) {
                            return candidatos;
                        }
                    }
                    
                    if(cnt >= T) {
                        return candidatos;
                    }
                }
            }
            return candidatos;
        }
        
        vector<tuple<int, vector<double>, double>> C_ANN_K(const vector<double>& query, double c, int k){
            //double t = 1.0;
            // if (N < 70000) t = 200.0;
            // else if (N < 500000) t = 1000.0;
            // else if (N < 2000000) t = 2000.0;
            // else t = 20000.0;
            // t *= 2.0;
            
            int T = 2*t*L + k;
            double r = R_min;
            
            vector<tuple<int, vector<double>, double>> acumulados;
            set<int> ids_usados;
            
            int rounds = 0;
            // const int MAX_ROUNDS = 30;
            
            while(true) {
                rounds++;
                auto nuevos = RC_NN_K(query, r, c, k, T);
                
                for(const auto& candidato : nuevos) {
                    int id = get<0>(candidato);
                    if(ids_usados.find(id) == ids_usados.end()) {
                        ids_usados.insert(id);
                        acumulados.push_back(candidato);
                    }
                }
                
                if((int)acumulados.size() >= k) {
                    sort(acumulados.begin(), acumulados.end(), 
                         [](const auto& a, const auto& b) { return get<2>(a) < get<2>(b); });
                    acumulados.resize(min(k, (int)acumulados.size()));
                    return acumulados;
                }
                
                r *= c;
            }
            
            sort(acumulados.begin(), acumulados.end(),
                 [](const auto& a, const auto& b) { return get<2>(a) < get<2>(b); });
            return acumulados;
        }

        int getDatasetSize() const { return datos.size(); }
};

vector<vector<double>> loadDataset(const string& path, size_t max_rows = 5000) {
    ifstream f(path);
    if (!f) throw runtime_error("No se pudo abrir " + path);

    string line;
    getline(f, line);
    vector<vector<double>> datos;
    datos.reserve(max_rows);

    while (getline(f, line) && datos.size() < max_rows) {
        stringstream ss(line);
        string cell;
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

// Métricas
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

double calcularRecallKNN(const vector<tuple<int, vector<double>, double>>& vecinos_dblsh,
                         const vector<pair<int, double>>& vecinos_reales) {
    if(vecinos_dblsh.empty() || vecinos_reales.empty()) return 0.0;
    
    set<int> ids_dblsh;
    for(const auto& v : vecinos_dblsh) {
        ids_dblsh.insert(get<0>(v));
    }
    
    set<int> ids_reales;
    for(const auto& v : vecinos_reales) {
        ids_reales.insert(v.first);
    }
    
    int interseccion = 0;
    for(int id : ids_dblsh) {
        if(ids_reales.count(id)) {
            interseccion++;
        }
    }
    
    int k = ids_reales.size();
    return (double)interseccion / k;
}

#define MAIN_C 1.5
#define MAIN_t 8000
#define MAIN_K 83
#define MAIN_L 2

int main(){
    std::filesystem::create_directories("results");

    cout << "============================================================" << endl;
    cout << "DB-LSH: Replicación de Gráficos del Paper (Fig. 5, 6, 7)" << endl;
    cout << "Variando n (proporción del dataset)" << endl;
    cout << "============================================================" << endl;

    const int D = 784;
    const double C = MAIN_C;
    const int t = MAIN_t;
    const int K = MAIN_K;
    const int L = MAIN_L;
    const int K_QUERIES = 50;  // Número de queries para promediar
    const int K_NN = 50;       // k vecinos a buscar
    const double R_MIN = 1;

    // Cargar dataset completo (60k)
    cout << "Cargando Fashion-MNIST completo...\n";
    vector<vector<double>> full_dataset = loadDataset("fashion_mnist.csv", 60000);
    cout << "Total cargado: " << full_dataset.size() << " imágenes" << endl;

    // Shuffle dataset
    mt19937 gen(42);
    shuffle(full_dataset.begin(), full_dataset.end(), gen);

    // Separar queries (primeras 50)
    vector<vector<double>> queries;
    for(int i = 0; i < K_QUERIES; i++) {
        queries.push_back(full_dataset[i]);
    }

    // Proporciones n a probar (como en Fig. 5, 6, 7)
    vector<double> n_values = {0.2, 0.4, 0.6, 0.8, 1.0};
    
    // Abrir archivo CSV para resultados
    ofstream csv_file("results/varying_n_results.csv");
    csv_file << "n,avg_recall,avg_ratio,avg_time_ms\n";
    
    cout << "\n" << string(70, '=') << endl;
    cout << "EXPERIMENTO: Variando n (proporción del dataset)" << endl;
    cout << string(70, '=') << endl;
    cout << "\nParámetros fijos:" << endl;
    cout << "  K (k-NN) = " << K_NN << endl;
    cout << "  Queries = " << K_QUERIES << endl;
    cout << "  C = " << C << ", R_min = " << R_MIN << ", t = " << t << endl;
    cout << "\nProbando n = {0.2, 0.4, 0.6, 0.8, 1.0}" << endl;
    cout << string(70, '-') << endl;
    
    for(double n : n_values) {
        // Calcular tamaño del dataset
        int dataset_size = static_cast<int>(n * (full_dataset.size() - K_QUERIES));
        
        cout << "\n[n = " << n << "] Dataset: " << dataset_size << " puntos" << endl;
        
        // Tomar subset del dataset (después de las queries)
        vector<vector<double>> dataset_index;
        for(int i = 0; i < dataset_size; i++) {
            dataset_index.push_back(full_dataset[K_QUERIES + i]);
        }
        
        // Construir índice DB-LSH
        cout << "  Construyendo índice DB-LSH..." << flush;
        DBLSH<K> indice(D, L, C, R_MIN, t, 42);
        indice.insertar(dataset_index);
        cout << " OK" << endl;
        
        // Ejecutar queries y medir métricas
        double total_recall = 0.0;
        double total_ratio = 0.0;
        double total_time_ms = 0.0;
        
        cout << "  Ejecutando " << K_QUERIES << " queries..." << flush;
        for(int q = 0; q < K_QUERIES; q++) {
            const auto& query = queries[q];
            
            // Medir tiempo de query
            auto start = high_resolution_clock::now();
            auto vecinos_dblsh = indice.C_ANN_K(query, C, K_NN);
            auto end = high_resolution_clock::now();
            
            double query_time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
            
            // Ground truth
            auto vecinos_reales = indice.encontrarKVecinosReales(query, K_NN);
            
            // Calcular métricas
            double recall = calcularRecallKNN(vecinos_dblsh, vecinos_reales);
            double ratio = calcularOverallRatioKNN(vecinos_dblsh, vecinos_reales);
            
            total_recall += recall;
            total_ratio += ratio;
            total_time_ms += query_time_ms;
        }
        cout << " OK" << endl;
        
        // Promedios
        double avg_recall = total_recall / K_QUERIES;
        double avg_ratio = total_ratio / K_QUERIES;
        double avg_time_ms = total_time_ms / K_QUERIES;
        
        cout << "  Resultados promedio:" << endl;
        cout << "    Recall: " << (avg_recall * 100) << "%" << endl;
        cout << "    Ratio: " << avg_ratio << "x" << endl;
        cout << "    Query time: " << avg_time_ms << " ms" << endl;
        
        // Guardar en CSV
        csv_file << n << "," << avg_recall << "," << avg_ratio << "," << avg_time_ms << "\n";
    }
    
    csv_file.close();
    
    cout << "\n" << string(70, '=') << endl;
    cout << "Resultados guardados en: results/varying_n_results.csv" << endl;
    cout << string(70, '=') << endl;
    
    // Mostrar tabla de resultados
    cout << "\nTabla de resultados:\n" << endl;
    cout << "   n    Recall(%)    Ratio    Time(ms)" << endl;
    cout << string(70, '-') << endl;
    
    ifstream csv_read("results/varying_n_results.csv");
    string line;
    getline(csv_read, line);  // Skip header
    
    while(getline(csv_read, line)) {
        stringstream ss(line);
        string n_str, recall_str, ratio_str, time_str;
        getline(ss, n_str, ',');
        getline(ss, recall_str, ',');
        getline(ss, ratio_str, ',');
        getline(ss, time_str, ',');
        
        double recall_pct = stod(recall_str) * 100;
        cout << "  " << n_str << "    " 
             << recall_pct << "%    " 
             << ratio_str << "    " 
             << time_str << endl;
    }
    csv_read.close();

    return 0;
}
