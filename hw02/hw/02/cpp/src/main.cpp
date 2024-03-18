#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Plane_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Plane_3 Plane_3;

using json = nlohmann::json;

void  generate_lod02(json &j);
void  generate_lod12(json &j);


int main(int argc, const char * argv[]) {
  //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
  const char* filename = (argc > 1) ? argv[1] : "../../data/twobuildings.city.json";
  std::cout << "Processing: " << filename << std::endl;
  std::ifstream input(filename);
  json j;
  input >> j; //-- store the content of the file in a nlohmann::json object
  input.close();

  generate_lod02(j);
  generate_lod12(j);
  //-- write to disk the modified city model (out.city.json)
  std::ofstream o("out.city.json");
  o << j.dump(2) << std::endl;
  o.close();

  return 0;
}

//// orient the surface by the area of the outer ring and the inner ring.
void orient_test(json& j, json& surface) {
    Polygon_2  outer_ring;
    auto& outer_boundary = surface[0];
    for (auto& vi : outer_boundary) {
        std::vector<int> v = j["vertices"][vi.get<int>()];
        double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
        double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
        outer_ring.push_back(K::Point_2(x, y));
    }
    if (outer_ring.area() > 0) {
        std::reverse(outer_boundary.begin(), outer_boundary.end());
    }

    if (surface.size() > 1) {
        for (int i = 1; i < surface.size(); i++) {
            auto& inner_boundary = surface[i];
            Polygon_2 inner_ring;
            for (auto& vi : inner_boundary) {
                std::vector<int> v = j["vertices"][vi.get<int>()];
                double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                inner_ring.push_back(K::Point_2(x, y));
            }
            if (inner_ring.area() < 0) {
                std::reverse(surface[i].begin(), surface[i].end());
            }
        }
    }
}

//// remove the boundary vertices with 10% of the highest z value.
//// find the highest z value/index of the left vertices as the roof point.
std::pair<int, double> roof_point(std::vector<std::pair<double, int>> height_index,double rate = 0.1){
    std::sort(height_index.begin(), height_index.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first; //
              });
    int remove_index = height_index.size() * rate; // remove 10% of the highest z value.
    if (remove_index > 0) {
        height_index.erase(height_index.begin(), height_index.begin() + remove_index);
    }
    //find the highest z value/index of the left vertices.
    double max_z = -INFINITY;
    int max_index = -1;
    for (auto &hi: height_index) {
        if (hi.first > max_z) {
            max_z = hi.first;
            max_index = hi.second;
        }
    }
    // std::cout << "max:" << max_z << std::endl;
    return std::make_pair(max_index, max_z);
}

//// record the height of the vertices and their index.
std::vector<std::pair<double, int>> get_height_index(json& first_geo,json& j){
    std::vector<std::pair<double, int>> height_index;
    if (first_geo["type"] == "Solid"){
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            for (int k = 0; k < first_geo["boundaries"][i].size(); k++) {
                // get the z value of the vertices and store them in the height vector.
                for (auto &vi_list: first_geo["boundaries"][i][k]) {
                    for (auto &vi: vi_list) {
                        int vertex_index = vi.get<int>();
                        auto &vertex = j["vertices"][vertex_index];
                        int vertex_z = vertex[2].get<int>();
                        double z = (vertex_z * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                        height_index.push_back(std::make_pair(z, vertex_index));
                    }
                }
            }
        }
    }
    else if (first_geo["type"] == "MultiSurface") {
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            for (auto &vi_list: first_geo["boundaries"][i]) {
                for (auto &vi: vi_list) {
                    int vertex_index = vi.get<int>();
                    auto &vertex = j["vertices"][vertex_index];
                    int vertex_z = vertex[2].get<int>();
                    double z = (vertex_z * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                    height_index.push_back(std::make_pair(z, vertex_index));
                }
            }
        }
    }
    return height_index;
}

//// select the ground surface using mean z-value, instead of referring to semantics.
//// ground surface is identified by the lowest mean z-value and the largest area.
void get_ground_surface(json& j, json& first_geo, json& ground_surface, json& ground_surface_sem_index) {
    int shell_index = -1;
    int ring_index = -1;
    double z_mean = INFINITY;
    double area = -INFINITY;

    if (first_geo["type"] == "Solid") {
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            for (int k = 0; k < first_geo["boundaries"][i].size(); k++) {
                double surface_sum_z = 0;
                int surface_num_vi = 0;
                Polygon_2 surface_2d; // contains x and y coordinates of the outer boundary of the original surface.
                auto& outer_boundary = first_geo["boundaries"][i][k][0];
                for (auto& vi : outer_boundary) {
                    std::vector<int> v = j["vertices"][vi.get<int>()];
                    double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                    double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                    double z = (v[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                    surface_2d.push_back(K::Point_2(x, y));
                    surface_sum_z += z;
                    surface_num_vi += 1;
                }
                double surface_mean_z = surface_sum_z / surface_num_vi;
                double surface_2d_area = fabs(surface_2d.area());
                if (surface_mean_z < z_mean && surface_2d_area >= area) {
                    shell_index = i;
                    ring_index = k;
                    z_mean = surface_mean_z;
                    area = surface_2d_area;
                }
            }
        }
        orient_test(j, first_geo["boundaries"][shell_index][ring_index]);
        ground_surface.push_back(first_geo["boundaries"][shell_index][ring_index]);
        ground_surface_sem_index.push_back(0);
    }
    else if (first_geo["type"] == "MultiSurface") {
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            double surface_sum_z = 0;
            int surface_num_vi = 0;
            Polygon_2 surface_2d;
            auto& outer_boundary = first_geo["boundaries"][i][0];
            for (auto& vi : outer_boundary) {
                std::vector<int> v = j["vertices"][vi.get<int>()];
                double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                double z = (v[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                surface_2d.push_back(K::Point_2(x, y));
                surface_sum_z += z;
                surface_num_vi += 1;
            }
            double surface_mean_z = surface_sum_z / surface_num_vi;
            double surface_2d_area = fabs(surface_2d.area());
            if (surface_mean_z <= z_mean && surface_2d_area >= area) {
                ring_index = i;
                z_mean = surface_mean_z;
                area = surface_2d_area;
            }
        }
        orient_test(j, first_geo["boundaries"][ring_index]);
        ground_surface.push_back(first_geo["boundaries"][ring_index]);
        ground_surface_sem_index.push_back(0);
    }
}

//// select the surface containing the roof point as the roof surfaces
//// and the flattest surface as the main roof surface.
auto get_roof_surface(json& first_geo, json& j, int roof_index){
    if (first_geo["type"] == "Solid"){
        std::vector<std::pair<int, int>> surface_index;
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            for (int k = 0; k < first_geo["boundaries"][i].size(); k++) {
                for (auto &vi_list: first_geo["boundaries"][i][k]) {
                    for (auto &vi: vi_list) {
                        if (vi.get<int>() == roof_index) {
                            surface_index.emplace_back(i, k);
                        }
                    }
                }
            }
        }
        // calculate the normal vector of the surface and select the flattest surface as the main roof surface.
        double roof_normal = -INFINITY;
        int roof_index = -1;
        for (int i = 0; i < surface_index.size(); i++) {
            int shell_index = surface_index[i].first;
            int ring_index = surface_index[i].second;
            std::vector<Point_3> surface_3d;
            for (auto &vi_list: first_geo["boundaries"][shell_index][ring_index]) {
                for (auto &vi: vi_list) {
                    std::vector<int> v = j["vertices"][vi.get<int>()];
                    double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                    double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                    double z = (v[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                    surface_3d.emplace_back(x, y, z);
                }
            }
            Vector_3 normal = CGAL::unit_normal(surface_3d[0], surface_3d[1], surface_3d[2]);
            double normal_z = fabs(normal.z()); // avoid negative value caused by wrong orientation.
            if (normal_z > roof_normal) {
                roof_normal = normal_z;
                roof_index = i;
            }
        }
        //std::cout << "roof_normal:" << roof_normal << std::endl;
        auto main_roof = first_geo["boundaries"][surface_index[roof_index].first][surface_index[roof_index].second];
        return main_roof;
    }
    else if (first_geo["type"] == "MultiSurface") {
        std::vector<int> surface_index;
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            for (auto &vi_list: first_geo["boundaries"][i]) {
                for (auto &vi: vi_list) {
                    if (vi.get<int>() == roof_index) {
                        surface_index.push_back(i);
                    }
                }
            }
        }

        double roof_normal = -INFINITY;
        int roof_index = -1;
        for (auto & i : surface_index) {
            int ring_index = i;
            std::vector<Point_3> surface_3d;
            for (auto &vi_list: first_geo["boundaries"][ring_index]) {
                for (auto &vi: vi_list) {
                    std::vector<int> v = j["vertices"][vi.get<int>()];
                    double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                    double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                    double z = (v[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                    surface_3d.push_back(Point_3(x, y, z));
                }
            }
            Vector_3 normal = CGAL::unit_normal(surface_3d[0], surface_3d[1], surface_3d[2]);
            double normal_z = normal.z();
            if (normal_z > roof_normal) {
                roof_normal = normal_z;
                roof_index = ring_index;
            }
        }
        auto main_roof = first_geo["boundaries"][roof_index];
        return main_roof;
    }
}

//// generate LoD0.2 geometry
void generate_lod02(json& j) {
    for (auto& co : j["CityObjects"].items()) {
        // if an CityObject already has lod0.2, do nothing and skip.
        bool lod02_exist = false;
        for (auto& g : co.value()["geometry"]) {
            if (g["lod"] == "0.2") {
                lod02_exist = true;
                break;
            }
        }
        if (lod02_exist) {
            std::cout << "This CityObject already contains LoD0.2 geometry, skip." << std::endl;
            continue;
        }

        if (!co.value()["geometry"].empty()) {
            // initialize ground surface array and ground surface semantic index array.
            auto ground_surface = json::array();
            auto ground_surface_sem_index = json::array();

            // get the first geometry which has the type "Solid" or "MultiSurface", use this to generate lod0.2.
            json first_geo;
            for (auto& g : co.value()["geometry"]) {
                if (g["type"] == "Solid" | g["type"] == "MultiSurface") {
                    first_geo = g;
                    break;
                }
            }

            // get ground surfaces using mean z-value, instead of referring to semantics.
            get_ground_surface(j, first_geo, ground_surface, ground_surface_sem_index);

            // create new lod0.2 geometry
            json new_geometry = {{"lod", "0.2"}, {"type", "MultiSurface"}};
            auto sem_array = json::array({{{"type", "GroundSurface"}}});
            new_geometry["semantics"]["surfaces"] = sem_array;
            new_geometry["semantics"]["values"] = ground_surface_sem_index;
            new_geometry["boundaries"] = ground_surface;
            co.value()["geometry"].push_back(new_geometry);
        }
    }
}

//// generate LoD1.2 geometry
void generate_lod12(json& j){
    for (auto& co : j["CityObjects"].items()) {
        // if an CityObject already has lod1.2, do nothing and skip.
        bool lod12_exist = false;
        for (auto &g: co.value()["geometry"]) {
            if (g["lod"] == "1.2") {
                lod12_exist = true;
                break;
            }
        }
        if (lod12_exist) {
            std::cout << "This CityObject already contains LoD0.2 geometry, skip." << std::endl;
            continue;
        }
        if (!co.value()["geometry"].empty()) {
            // surface12/index12:: the exterior shell of each object.
            auto surface12 = json::array();
            auto surface_index12 = json::array();
            // ground_surface12/index12:: the ground surface of each object.
            auto shell_surface12 = json::array();
            auto shell_surface_index12 = json::array();
            auto ground_surface12 = json::array();
            // get the ground surface from lod0.2 geometry.
            for (auto &g:co.value()["geometry"]) {
                if (g["lod"] == "0.2") {
                    ground_surface12 = g["boundaries"]; // used for iterating later
                    shell_surface12 = g["boundaries"];
                    shell_surface_index12 = g["semantics"]["values"];
                    break;
                }
            }
            json first_geo;
            // generate the roof surface
            for (auto &g: co.value()["geometry"]) {
                if (g["type"] == "Solid" | g["type"] == "MultiSurface") {
                    first_geo = g;
                    break;
                }
            }
            // record the height of the vertices and their index.
            std::vector<std::pair<double, int>> height_index = get_height_index(first_geo, j);
            int roof_index = roof_point(height_index,0.1).first;
            double roof_z = roof_point(height_index,0.1).second;
            auto main_roof = get_roof_surface(first_geo, j, roof_index);
            // get the lowest point of the roof surface as the eave height.
            double eave_z = INFINITY;
            for (auto &vi_list: main_roof) {
                for (auto &vi: vi_list) {
                    int vertex_index = vi.get<int>();
                    auto &vertex = j["vertices"][vertex_index];
                    int vertex_z = vertex[2].get<int>();
                    double z = (vertex_z * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                    if (z < eave_z) {
                        eave_z = z;
                    }
                }
            }
            double roof_height12 = eave_z + 0.7*(roof_z - eave_z);


            //// create new lod1.2 geometry
            json new_geometry12 = {{"lod", "1.2"}, {"type", "Solid"}};
            auto sem_array12 = json::array({{{"type", "GroundSurface"}},
                                            {{"type", "RoofSurface"}},
                                            {{"type", "WallSurface"}}});
            new_geometry12["semantics"]["surfaces"] = sem_array12;

            // use ground surface to generate roof surface
            for (auto& gs: ground_surface12){
                auto roof_surface = json::array();
                auto& vertices = j["vertices"];
                auto& transform = j["transform"];

                for (auto& index_list : gs) {
                    auto rf_lift = json::array();
                    for (auto& index: index_list) {
                        int vertex_index = index.get<int>();
                        std::vector<int> vi = vertices[vertex_index];
                        double z_lift = roof_height12;
                        int z_int = std::round((z_lift - transform["translate"][2].get<double>()) / transform["scale"][2].get<double>());
                        vertices.push_back(json::array({vi[0], vi[1], z_int}));
                        rf_lift.push_back(vertices.size() - 1);
                    }
                    // roof orientation is reversed of the ground surface.
                    std::reverse(rf_lift.begin(), rf_lift.end());
                    roof_surface.push_back(rf_lift);
                    shell_surface12.push_back(roof_surface);
                    shell_surface_index12.push_back(1);
                }
                surface12.push_back(shell_surface12);
                surface_index12.push_back(shell_surface_index12);
            }
            new_geometry12["boundaries"] = surface12;
            new_geometry12["semantics"]["values"] = surface_index12;
            co.value()["geometry"].push_back(new_geometry12);
        }
    }
}