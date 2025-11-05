#include "R_star.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

// Constructor
RStarTreeIndex::RStarTreeIndex() : total_elements_(0) {}

void RStarTreeIndex::insertPrueba(int id , tuple<double,double> data) {
    Point p;
    bg::set<0>(p, get<0>(data));
    bg::set<1>(p, get<1>(data));
    rtree_.insert(make_pair(p, id));
    total_elements_++;
}


bool RStarTreeIndex::loadPrueba(vector<tuple<double,double>>  filepath) {
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

// Window Query: buscar puntos dentro de un rectángulo
vector<RStarTreeIndex::Value> RStarTreeIndex::windowQuery(double x_min, double y_min, double x_max, double y_max) const {
    Box query_box(Point(x_min, y_min), Point(x_max, y_max));
    
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

