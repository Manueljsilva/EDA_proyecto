#include "R_star.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

// Constructor
RStarTreeIndex::RStarTreeIndex() : total_elements_(0) {}

void RStarTreeIndex::insertPrueba(int id , array<double,10> data) {
    Point p;
    bg::set<0>(p, data[0]);
    bg::set<1>(p, data[1]);
    bg::set<2>(p, data[2]);
    bg::set<3>(p, data[3]);
    bg::set<4>(p, data[4]);
    bg::set<5>(p, data[5]);
    bg::set<6>(p, data[6]);
    bg::set<7>(p, data[7]);
    bg::set<8>(p, data[8]);
    bg::set<9>(p, data[9]);
    rtree_.insert(make_pair(p, id));
    total_elements_++;
}


bool RStarTreeIndex::loadPrueba(vector<array<double,10>>  filepath) {
    int id = 1 ;
    for (auto linea : filepath) {    
        insertPrueba(id , linea);
        id++;
    }
    cout << "Cargados " << id-1 << " registros al R*-tree" << endl;
    return true;
}
// Imprimir estadísticas
void RStarTreeIndex::printStats() const {
    cout << "  Estadísticas del R*-tree:" << endl;
    cout << "   Total de elementos: " << total_elements_ << endl;
    cout << "   Parámetro R*: 16 (max elementos por nodo)" << endl;
}

// Limpiar
void RStarTreeIndex::clear() {
    rtree_.clear();
    total_elements_ = 0;
}

// Window Query: buscar puntos dentro de un hiper-rectángulo 10D
vector<RStarTreeIndex::Value> RStarTreeIndex::windowQuery(const array<double,10>& mins, const array<double,10>& maxs) const {
    Point p_min, p_max;
    
    // Construir puntos min y max usando los arrays
    bg::set<0>(p_min, mins[0]); bg::set<0>(p_max, maxs[0]);
    bg::set<1>(p_min, mins[1]); bg::set<1>(p_max, maxs[1]);
    bg::set<2>(p_min, mins[2]); bg::set<2>(p_max, maxs[2]);
    bg::set<3>(p_min, mins[3]); bg::set<3>(p_max, maxs[3]);
    bg::set<4>(p_min, mins[4]); bg::set<4>(p_max, maxs[4]);
    bg::set<5>(p_min, mins[5]); bg::set<5>(p_max, maxs[5]);
    bg::set<6>(p_min, mins[6]); bg::set<6>(p_max, maxs[6]);
    bg::set<7>(p_min, mins[7]); bg::set<7>(p_max, maxs[7]);
    bg::set<8>(p_min, mins[8]); bg::set<8>(p_max, maxs[8]);
    bg::set<9>(p_min, mins[9]); bg::set<9>(p_max, maxs[9]);
    
    Box query_box(p_min, p_max);
    
    vector<Value> result;
    rtree_.query(bgi::intersects(query_box), back_inserter(result));
    
    return result;
}

// Window Query alternativo: usando un Box directamente
vector<RStarTreeIndex::Value> RStarTreeIndex::windowQuery(const Box& query_box) const {
    vector<Value> result;
    rtree_.query(bgi::intersects(query_box), back_inserter(result));
    
    return result;
}

