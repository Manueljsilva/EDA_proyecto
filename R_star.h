#ifndef RSTARTREEINDEX_H
#define RSTARTREEINDEX_H

#include <vector>
#include <string>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class RStarTreeIndex {
public:
    // CAMBIAR EN CASO DE OTRO DATASET
    using Point = bg::model::point<double, 2, bg::cs::cartesian>;
    using Box = bg::model::box<Point>;
    using Value = pair<Point, int>; // par de punto y ID
    
private:
    // R*-tree con parámetro 16 (máximo de elementos por nodo)
    bgi::rtree<Value, bgi::rstar<16>> rtree_;
    int total_elements_;
    
public:
    // Constructor
    RStarTreeIndex();

    void insertPrueba(int id , tuple<double,double> data);
    bool loadPrueba(vector<tuple<double,double>> filepath);

    // Window Query: buscar puntos dentro de un rectángulo
    vector<Value> windowQuery(double x_min, double y_min, double x_max, double y_max) const;
    
    // Window Query alternativo: usando un Box directamente
    vector<Value> windowQuery(const Box& query_box) const;

    // Limpiar el índice
    void clear();
    // Imprimir estadísticas
    void printStats() const;
private:
};

#endif // RSTARTREEINDEX_H
