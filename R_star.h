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
    using Point = bg::model::point<double, 10, bg::cs::cartesian>;
    using Box = bg::model::box<Point>;
    using Value = pair<Point, int>; // par de punto y ID
    
private:
    // R*-tree con parámetro 16 (máximo de elementos por nodo)
    bgi::rtree<Value, bgi::rstar<16>> rtree_;
    int total_elements_;
    
public:
    // Constructor
    RStarTreeIndex();

    void insertPrueba(int id , array<double,10> data);
    bool loadPrueba(vector<array<double,10>> filepath);

    // Window Query: buscar puntos dentro de un hiper-rectángulo 10D
    vector<Value> windowQuery(const array<double,10>& mins, const array<double,10>& maxs) const;
    
    // Window Query alternativo: usando un Box directamente
    vector<Value> windowQuery(const Box& query_box) const;

    // Limpiar el índice
    void clear();
    // Imprimir estadísticas
    void printStats() const;
private:
};

#endif // RSTARTREEINDEX_H
