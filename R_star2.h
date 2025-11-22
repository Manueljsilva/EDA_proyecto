#ifndef RSTARTREEINDEX_H
#define RSTARTREEINDEX_H

#include <vector>
#include <string>
#include <array>
#include <type_traits>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template <size_t Dim>
class RStarTreeIndex {
public:
    using Point = bg::model::point<double, Dim, bg::cs::cartesian>;
    using Box = bg::model::box<Point>;
    using Value = pair<Point, int>; // par de punto y ID
    
private:
    // R*-tree con parámetro 16 (máximo de elementos por nodo)
    bgi::rtree<Value, bgi::rstar<16>> rtree_;
    int total_elements_;

    template <size_t I = 0>
    static std::enable_if_t<I == Dim> fillPoint(Point&, const std::array<double, Dim>&) {}

    template <size_t I = 0>
    static std::enable_if_t<I < Dim> fillPoint(Point& p, const std::array<double, Dim>& data) {
        bg::set<I>(p, data[I]);
        fillPoint<I + 1>(p, data);
    }
    
public:
    // Constructor
    RStarTreeIndex() : total_elements_(0) {};

    void insertPrueba(int id , array<double,Dim> data) {
        Point p;
        fillPoint(p, data);
        rtree_.insert(std::make_pair(p, id));
        total_elements_++;
    }
    
    // Bulk-loading: construccion eficiente del R*-tree (Dim-dimensional)
    void bulkLoad(const vector<pair<array<double, Dim>, int>>& data) {
        vector<Value> values;
        values.reserve(data.size());

        for (const auto& item : data) {
            const auto& hash_array = item.first;
            int id = item.second;
            Point p;
            fillPoint(p, hash_array);
            values.push_back(std::make_pair(p, id));
        }

        // Reconstruir R*-tree usando bulk-loading (mas eficiente)
        rtree_.clear();
        rtree_ = bgi::rtree<Value, bgi::rstar<16>>(values.begin(), values.end());
        total_elements_ = static_cast<int>(values.size());
    }
    
    bool loadPrueba(vector<array<double,Dim>> filepath) {
        int id = 1 ;
        for (auto linea : filepath) {    
            insertPrueba(id , linea);
            id++;
        }
        cout << "Cargados " << id-1 << " registros al R*-tree" << endl;
        return true;
    }

    // Window Query: buscar puntos dentro de un hiper-rectángulo
    vector<Value> windowQuery(const array<double,Dim>& mins, const array<double,Dim>& maxs) const {
        Point p_min, p_max;
        fillPoint(p_min, mins);
        fillPoint(p_max, maxs);
        Box query_box(p_min, p_max);
        std::vector<Value> result;
        rtree_.query(bgi::intersects(query_box), std::back_inserter(result));
        return result;
    }
    
    // Window Query alternativo: usando un Box directamente
    vector<Value> windowQuery(const Box& query_box) const {
        vector<Value> result;
        rtree_.query(bgi::intersects(query_box), back_inserter(result));
        
        return result;
    }

    // Limpiar el índice
    void clear() {
        rtree_.clear();
        total_elements_ = 0;
    }

    // Imprimir estadísticas
    void printStats() const {
        cout << "  Estadísticas del R*-tree:" << endl;
        cout << "   Total de elementos: " << total_elements_ << endl;
        cout << "   Dimensiones: " << Dim << "D" << endl;
        cout << "   Parámetro R*: 16 (max elementos por nodo)" << endl;
    }
};

#endif // RSTARTREEINDEX_H
