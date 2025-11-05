#include <iostream>
#include <fstream>
#include "R_star.h"
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
using namespace std;


vector <vector <double>> a = {
    {0.6, 0.8},
    {0.3, -0.9}
};
// clase  
class DBfsh{
    private:
        int D = 2; // dimensiones
        int L = 5; // numero de tablas hash
        int K = 12; // numero de funciones hash 
        double C = 1.5; // const de error
        double w0 = 0; // ancho de la ventana
        int t = 1; // parámetro t (entero positivo)
        RStarTreeIndex indice;
        
        // Mapeo: hash -> punto original (para recuperar datos originales)
        vector<pair<tuple<double,double>, tuple<double,double>>> hash_to_original;


        vector<tuple<double,double>> funcionHash(vector<tuple<double,double>> punto){
            // implementación de la función hash
            vector <tuple<double,double>> proyecciones;
            for(auto[x,y]: punto){
                double h1 = a[0][0]*x + a[0][1]*y;
                double h2 = a[1][0]*x + a[1][1]*y;
                proyecciones.push_back({h1,h2});
            }
            return proyecciones;
        }
        tuple<double,double> funcionHash(tuple<double,double> punto){
            // implementación de la función hash
            double x = get<0>(punto);
            double y = get<1>(punto);
            double h1 = a[0][0]*x + a[0][1]*y;
            double h2 = a[1][0]*x + a[1][1]*y;
            return {h1,h2};
        }
        // distancia euclidiana
        double distanciaEuclidiana(tuple<double,double> p1, tuple<double,double> p2){
            double dx = get<0>(p1) - get<0>(p2);
            double dy = get<1>(p1) - get<1>(p2);
            return sqrt(dx*dx + dy*dy);
        }


    public:
        // constructor por defecto
        DBfsh(){
            w0 = 4*C*C;
        }
        // constructor con parámetros
        DBfsh(int dim, int L_, int K_, double C_ , int t_) : D(dim), L(L_), K(K_), C(C_), w0(4*C*C), t(t_) {}

        void insertar(vector<tuple<double , double>> valor){
            vector<tuple<double,double>> PuntosHash; 
            hash_to_original.clear(); // Limpiar mapeo anterior
            
            for (auto punto : valor){
                tuple<double,double> hash_punto = funcionHash(punto);
                PuntosHash.push_back(hash_punto);
                // Guardar mapeo: hash -> original
                hash_to_original.push_back({hash_punto, punto});
            }
            
            cout << "Hashes generados:" << endl;
            for (size_t i = 0; i < PuntosHash.size(); i++) {
                cout << "Original: (" << get<0>(valor[i]) << ", " << get<1>(valor[i]) << ") -> "
                     << "Hash: (" << get<0>(PuntosHash[i]) << ", " << get<1>(PuntosHash[i]) << ")\n";
            }
            indice.loadPrueba(PuntosHash);
            cout << endl;
        }
        void imprimir(){
            indice.printStats();
        }
        
        
        // Algorithm 1: (r,c)-NN Query
        // Input: q (query point), r (query radius), c (approximation ratio), t (positive integer)
        // Output: A point o or ∅
        tuple<double,double> RC_NN(tuple<double,double> query, double r, double c){
            int cnt = 0;
            
            cout << "\n--- (r,c)-NN Query ---" << endl;
            cout << "Query: (" << get<0>(query) << ", " << get<1>(query) << ")" << endl;
            cout << "r = " << r << ", c = " << c << ", t = " << t << ", L = " << L << endl;
            
            // Calcular G_i(q) para cada tabla (en nuestro caso simplificado, solo 1 tabla)
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
                    tuple<double,double> hash_candidato = {bg::get<0>(res.first), bg::get<1>(res.first)};
                    
                    // Buscar el punto original correspondiente a este hash
                    tuple<double,double> punto_original = {0.0, 0.0};
                    for(const auto& [hash, original] : hash_to_original) {
                        if(abs(get<0>(hash) - get<0>(hash_candidato)) < 0.0001 && 
                           abs(get<1>(hash) - get<1>(hash_candidato)) < 0.0001) {
                            punto_original = original;
                            break;
                        }
                    }
                    
                    cnt = cnt + 1;
                    cout << "  Punto " << cnt << ": hash=(" << get<0>(hash_candidato) << ", " << get<1>(hash_candidato) << ")";
                    cout << " -> original=(" << get<0>(punto_original) << ", " << get<1>(punto_original) << ")" << endl;
                    
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
            return {INFINITY, INFINITY}; // Representamos ∅ con infinito
        }
        
        // Algorithm 2: c-ANN Query
        // Input: q (query point), c (approximation ratio)
        // Output: A point o
        tuple<double,double> C_ANN(tuple<double,double> query, double c){
            cout << "\n" << string(60, '=') << endl;
            cout << "c-ANN Query" << endl;
            cout << string(60, '=') << endl;
            cout << "Query: (" << get<0>(query) << ", " << get<1>(query) << ")" << endl;
            cout << "c = " << c << endl;
            
            double r = 1.0; // r ← 1
            int iteracion = 0;
            
            // while TRUE do
            while(true){
                iteracion++;
                cout << "\n--- Iteración " << iteracion << " ---" << endl;
                cout << "r = " << r << endl;
                
                // o ← call (r,c)-NN
                tuple<double,double> o = RC_NN(query, r, c);
                
                // if o ≠ ∅ then return o
                if(get<0>(o) != INFINITY) {
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
                    return {INFINITY, INFINITY};
                }
            }
        }


};



int main(){
    // dataset de dos dimenciones
    vector<tuple<double, double>> datos = {
        {1.0, 1.0},
        {2.0, 2.0},
        {4.0, 2.0},
        {5.0, 5.0},
        {7.0, 8.0}
    };

    cout << "=====================================" << endl;
    cout << "DB-LSH con (r,c)-NN y c-ANN Query" << endl;
    cout << "=====================================" << endl;
    cout << "Parámetros:" << endl;
    cout << "  L = 1 (tablas hash)" << endl;
    cout << "  K = 2 (funciones hash)" << endl;
    cout << "  C = 1.5 (constante de error)" << endl;
    cout << "  t = 1 (parámetro positivo)" << endl;
    cout << "  w0 = 4*C^2 = 9.0" << endl;
    cout << "=====================================" << endl;

    DBfsh indice(2, 1, 2, 1.5, 1);
    indice.insertar(datos);
    indice.imprimir();

    // Probar c-ANN Query (Algorithm 2)
    cout << "\n" << string(60, '=') << endl;
    cout << "PRUEBA DE c-ANN QUERY (Algorithm 2)" << endl;
    cout << string(60, '=') << endl;
    
    tuple<double, double> query = {6.0, 6.0};
    double c = 1.5;
    
    tuple<double,double> vecino = indice.C_ANN(query, c);
    
    cout << "\n" << string(60, '=') << endl;
    cout << "RESULTADO FINAL:" << endl;
    if(get<0>(vecino) != INFINITY) {
        cout << "  Vecino más cercano (c-aproximado): (" 
             << get<0>(vecino) << ", " << get<1>(vecino) << ")" << endl;
        
        // Calcular distancia real
        double dx = get<0>(query) - get<0>(vecino);
        double dy = get<1>(query) - get<1>(vecino);
        double dist = sqrt(dx*dx + dy*dy);
        cout << "  Distancia real: " << dist << endl;
    } else {
        cout << "  No se encontró vecino" << endl;
    }
    cout << string(60, '=') << endl;

    return 0;
}


    

    


