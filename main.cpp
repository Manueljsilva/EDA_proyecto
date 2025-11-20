#include <iostream>
#include <fstream>
#include "R_star.h"
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;
// clase  
class DBfsh{
    private:
        int D;      // Dimensión original (ej: 2, 10, 128)
        int K;      // Número de funciones hash = dimensión proyectada (siempre 2 para R*-tree 2D)
        int L;      // Número de tablas hash
        double C;   // Constante de aproximación
        double w0;  // Ancho de ventana base
        int t;      // Parámetro t (entero positivo)
        unsigned seed; // Semilla para reproducibilidad
        
        RStarTreeIndex indice;
        
        // Almacena los puntos originales (índice del vector = id del punto)
        vector<vector<double>> datos;

        // Matriz de proyección: K filas × D columnas
        // a[i][j] = coeficiente de la función hash i para la dimensión j
        vector<vector<double>> a;

        // Generar funciones hash aleatorias normalizadas
        void generarFuncionesHash() {
            a.resize(K);
            mt19937 gen(seed);
            normal_distribution<double> dist(0.0, 1.0);
            
            for(int i = 0; i < K; i++) {
                a[i].resize(D);
                double norm = 0.0;
                
                // Generar vector aleatorio N(0,1)
                for(int j = 0; j < D; j++) {
                    a[i][j] = dist(gen);
                    norm += a[i][j] * a[i][j];
                }
                
                // Normalizar ||a[i]|| = 1
                norm = sqrt(norm);
                for(int j = 0; j < D; j++) {
                    a[i][j] /= norm;
                }
            }
        }
        
        // Proyectar punto N-dimensional → K-dimensional (siempre K=2 para R*-tree)
        tuple<double,double> funcionHash(const vector<double>& punto) {
            if(punto.size() != D) {
                throw runtime_error("Punto debe tener " + to_string(D) + " dimensiones");
            }
            
            // h_i(p) = a[i] · p (producto punto)
            double h1 = 0.0, h2 = 0.0;
            
            for(int j = 0; j < D; j++) {
                h1 += a[0][j] * punto[j];
                h2 += a[1][j] * punto[j];
            }
            
            return {h1, h2};
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
        // Constructor con parámetros (K siempre = 2 para R*-tree 2D)
        DBfsh(int dim, int L_, double C_, int t_, unsigned seed_ = 42) 
            : D(dim), K(2), L(L_), C(C_), w0(4*C_*C_), t(t_), seed(seed_) {
            generarFuncionesHash();
            
            cout << "DB-LSH inicializado:" << endl;
            cout << "  Dimensión original: " << D << "D" << endl;
            cout << "  Dimensión proyectada: " << K << "D (fijo para R*-tree)" << endl;
            cout << "  Tablas hash: " << L << endl;
            cout << "  C = " << C << ", w0 = " << w0 << ", t = " << t << endl;
            cout << "  Semilla: " << seed << endl;
        }

        void insertar(const vector<vector<double>>& datos_input){
            // Guardar datos originales
            datos = datos_input;
            
            cout << "\nInsertando " << datos.size() << " puntos de " << D << "D..." << endl;
            
            // Proyectar y guardar en R*-tree con ID = índice del vector
            for (size_t i = 0; i < datos.size(); i++){
                tuple<double,double> hash_punto = funcionHash(datos[i]);
                int id = static_cast<int>(i);  // ID = índice en vector datos
                indice.insertPrueba(id, hash_punto);
            }
            
            cout << "Proyecciones generadas (primeros 5):" << endl;
            for (size_t i = 0; i < min((size_t)5, datos.size()); i++) {
                auto hash = funcionHash(datos[i]);
                cout << "  Punto[" << i << "] " << D << "D -> Hash: (" 
                     << get<0>(hash) << ", " << get<1>(hash) << ")" << endl;
            }
            if(datos.size() > 5) {
                cout << "  ... (" << (datos.size() - 5) << " más)" << endl;
            }
            cout << endl;
        }
        void imprimir(){
            indice.printStats();
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
            
            // Calcular G_i(q) para cada tabla
            tuple<double,double> hash_query = funcionHash(query);
            cout << "Hash G(q) = (" << get<0>(hash_query) << ", " << get<1>(hash_query) << ")" << endl;
            
            // Para cada tabla i = 1 to L
            for(int i = 0; i < L; i++){
                cout << "\nTabla " << (i+1) << ":" << endl;
                
                // while a point o ∈ W(G_i(q), w_0 · r) is found
                double w_r = w0 * r;
                double threshold = w_r / 2.0;
                
                double x_min = get<0>(hash_query) - threshold;
                double x_max = get<0>(hash_query) + threshold;
                double y_min = get<1>(hash_query) - threshold;
                double y_max = get<1>(hash_query) + threshold;
                
                cout << "  Ventana W(G(q), w_0·r = " << w_r << "): ";
                cout << "[" << x_min << ", " << x_max << "] x [" << y_min << ", " << y_max << "]" << endl;
                
                // Window Query en el R*-tree
                vector<RStarTreeIndex::Value> resultados = indice.windowQuery(x_min, y_min, x_max, y_max);
                cout << "  Puntos encontrados en ventana: " << resultados.size() << endl;
                
                // Procesar cada punto encontrado
                for(const auto& res : resultados) {
                    // El R*-tree ya nos da el ID directamente en res.second
                    int id = res.second;
                    
                    // Recuperar punto original usando el ID (O(1))
                    const vector<double>& punto_original = datos[id];
                    
                    cnt = cnt + 1;
                    cout << "  Punto " << cnt << ": id=" << id << endl;
                    
                    // if cnt = 2tL + 1 or ||q, o|| ≤ cr then return o;
                    double dist = distanciaEuclidiana(query, punto_original);
                    cout << "    dist(q, o) = " << dist << ", c·r = " << c*r << endl;
                    
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
            
            // return ∅
            cout << "  ∅ No se encontró punto satisfactorio" << endl;
            return vector<double>(D, INFINITY); // Representamos ∅ con infinitos
        }
        
        // Algorithm 2: c-ANN Query
        // Input: q (query point), c (approximation ratio)
        // Output: A point o
        vector<double> C_ANN(const vector<double>& query, double c){
            cout << "\n" << string(60, '=') << endl;
            cout << "c-ANN Query" << endl;
            cout << string(60, '=') << endl;
            cout << "Query " << D << "D: [";
            for(size_t i = 0; i < min((size_t)5, query.size()); i++) {
                cout << query[i];
                if(i < min((size_t)4, query.size()-1)) cout << ", ";
            }
            if(query.size() > 5) cout << " ...";
            cout << "]" << endl;
            cout << "c = " << c << endl;
            
            double r = 1.0; // r ← 1
            int iteracion = 0;
            
            // while TRUE do
            while(true){
                iteracion++;
                cout << "\n--- Iteración " << iteracion << " ---" << endl;
                cout << "r = " << r << endl;
                
                // o ← call (r,c)-NN
                vector<double> o = RC_NN(query, r, c);
                
                // if o ≠ ∅ then return o
                if(o[0] != INFINITY) {
                    cout << "\n✓ Punto encontrado!" << endl;
                    return o;
                }
                // else r ← cr
                else {
                    r = c * r;
                    cout << "  No encontrado, expandiendo: r ← c·r = " << r << endl;
                }
                
                // Límite de seguridad
                if(r > 1000.0) {
                    cout << "  ⚠ Alcanzado límite de expansión (r > 1000)" << endl;
                    return vector<double>(D, INFINITY);
                }
            }
        }


};



int main(){
    cout << "============================================================" << endl;
    cout << "DB-LSH Dinámico: N-Dimensional → 2D (R*-tree)" << endl;
    cout << "============================================================" << endl;
    
    // ============= EJEMPLO 1: 2D → 2D (caso base) =============
    cout << "\n### EJEMPLO 1: Datos 2D → Hash 2D ###\n" << endl;
    
    vector<vector<double>> datos_2d = {
        {1.0, 1.0},
        {2.0, 2.0},
        {4.0, 2.0},
        {5.0, 5.0},
        {7.0, 8.0}
    };

    DBfsh indice_2d(2, 1, 1.5, 1, 42);  // dim=2, L=1, C=1.5, t=1, seed=42
    indice_2d.insertar(datos_2d);
    indice_2d.imprimir();

    vector<double> query_2d = {6.0, 6.0};
    vector<double> vecino_2d = indice_2d.C_ANN(query_2d, 1.5);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO:" << endl;
    if(vecino_2d[0] != INFINITY) {
        cout << "  Vecino encontrado: [";
        for(size_t i = 0; i < vecino_2d.size(); i++) {
            cout << vecino_2d[i];
            if(i < vecino_2d.size()-1) cout << ", ";
        }
        cout << "]" << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;
    
    
    // ============= EJEMPLO 2: 5D → 2D =============
    cout << "\n\n### EJEMPLO 2: Datos 5D → Hash 2D ###\n" << endl;
    
    vector<vector<double>> datos_5d = {
        {1.0, 2.0, 3.0, 4.0, 5.0},
        {2.0, 3.0, 4.0, 5.0, 6.0},
        {5.0, 5.0, 5.0, 5.0, 5.0},
        {10.0, 10.0, 10.0, 10.0, 10.0},
        {1.5, 2.5, 3.5, 4.5, 5.5}
    };
    
    DBfsh indice_5d(5, 1, 1.5, 1, 123);  // dim=5, L=1, C=1.5, t=1, seed=123
    indice_5d.insertar(datos_5d);
    indice_5d.imprimir();
    
    vector<double> query_5d = {1.2, 2.1, 3.2, 4.1, 5.1};
    cout << "\nBuscando vecino cercano a query 5D..." << endl;
    vector<double> vecino_5d = indice_5d.C_ANN(query_5d, 1.5);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO:" << endl;
    if(vecino_5d[0] != INFINITY) {
        cout << "  Vecino encontrado 5D: [";
        for(size_t i = 0; i < vecino_5d.size(); i++) {
            cout << vecino_5d[i];
            if(i < vecino_5d.size()-1) cout << ", ";
        }
        cout << "]" << endl;
        
        // Calcular distancia real
        double dist = 0.0;
        for(size_t i = 0; i < query_5d.size(); i++) {
            double diff = query_5d[i] - vecino_5d[i];
            dist += diff * diff;
        }
        dist = sqrt(dist);
        cout << "  Distancia euclidiana: " << dist << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;
    
    
    // ============= EJEMPLO 3: 10D → 2D =============
    cout << "\n\n### EJEMPLO 3: Datos 10D → Hash 2D ###\n" << endl;
    
    vector<vector<double>> datos_10d;
    // Generar 10 puntos aleatorios en 10D
    for(int i = 0; i < 10; i++) {
        vector<double> punto(10);
        for(int j = 0; j < 10; j++) {
            punto[j] = (i + 1) * (j + 1) * 0.5;
        }
        datos_10d.push_back(punto);
    }
    
    DBfsh indice_10d(10, 1, 1.5, 1, 999);  // dim=10, L=1, C=1.5, t=1
    indice_10d.insertar(datos_10d);
    indice_10d.imprimir();
    
    vector<double> query_10d = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    cout << "\nBuscando vecino cercano en espacio 10D..." << endl;
    vector<double> vecino_10d = indice_10d.C_ANN(query_10d, 2.0);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO:" << endl;
    if(vecino_10d[0] != INFINITY) {
        cout << "  Vecino encontrado 10D (primeras 5 dims): [";
        for(size_t i = 0; i < min((size_t)5, vecino_10d.size()); i++) {
            cout << vecino_10d[i];
            if(i < 4) cout << ", ";
        }
        cout << " ...]" << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;
    
    cout << "\n✨ Demostración completa: LSH funciona con N dimensiones ✨" << endl;

    return 0;
}