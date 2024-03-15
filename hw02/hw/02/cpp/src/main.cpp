#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;

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


//// remove the boundary vertices with 10% of the highest z value.
int roof_point(std::vector<std::pair<double, int>> height_index){
    std::sort(height_index.begin(), height_index.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first; //
              });
    int remove_index = height_index.size() * 0.1; // remove 10% of the highest z value.
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
            return max_index;
        }
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

      // get the first geometry which has the type "Solid" or "MultiSurface",
      // use this to generate lod0.2.
      json first_geo;
      for (auto& g : co.value()["geometry"]) {
        if (g["type"] == "Solid" | g["type"] == "MultiSurface") {
          first_geo = g;
          break;
        }
      }

      // get ground surfaces using mean z-value, instead of referring to semantics.
      // ground surface is identified by the lowest mean z-value and the largest area.
      int shell_index = -1;
      int ring_index = -1;
      double z_mean = INFINITY;
      double area = -INFINITY;

      if (first_geo["type"] == "Solid") {
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
          for (int k = 0; k < first_geo["boundaries"][i].size(); k++) {
            double surface_sum_z = 0;
            int surface_num_vi = 0;
            Polygon_2 surface_2d; // this 2d surface only contains x and y coordinates of the original surface.

            for (auto&vi_list : first_geo["boundaries"][i][k]) {
              for (auto& vi : vi_list) {
                std::vector<int> v = j["vertices"][vi.get<int>()];
                double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                double z = (v[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                surface_2d.push_back(K::Point_2(x, y));
                surface_sum_z += z;
                surface_num_vi += 1;
              }
            }
            double surface_mean_z = surface_sum_z / surface_num_vi;
            double surface_2d_area = fabs(surface_2d.area());
            if (surface_mean_z < z_mean && surface_2d_area > area) {
              shell_index = i;
              ring_index = k;
              z_mean = surface_mean_z;
              area = surface_2d_area;
            }
          }
        }
      ground_surface.push_back(first_geo["boundaries"][shell_index][ring_index]);
      ground_surface_sem_index.push_back(0);
      }
      else if (first_geo["type"] == "MultiSurface") {
          for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            double surface_sum_z = 0;
            int surface_num_vi = 0;
            Polygon_2 surface_2d;

            for (auto&vi_list : first_geo["boundaries"][i]) {
              for (auto& vi : vi_list) {
                  std::vector<int> v = j["vertices"][vi.get<int>()];
                  double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                  double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                  double z = (v[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                  surface_2d.push_back(K::Point_2(x, y));
                  surface_sum_z += z;
                  surface_num_vi += 1;
              }
            }
            double surface_mean_z = surface_sum_z / surface_num_vi;
            double surface_2d_area = fabs(surface_2d.area());
            if (surface_mean_z <= z_mean && surface_2d_area >= area) {
              ring_index = i;
              z_mean = surface_mean_z;
              area = surface_2d_area;
            }
          }
          ground_surface.push_back(first_geo["boundaries"][ring_index]);
          ground_surface_sem_index.push_back(0);
      }

      // create new lod0.2 geometry
      json new_geometry = {{"lod", "0.2"}, {"type", "MultiSurface"}};
      auto sem_array = json::array({{{"type", "GroundSurface"}}, {{"type", "RoofSurface"}}});
      new_geometry["semantics"]["surfaces"] = sem_array;
      new_geometry["semantics"]["values"] = ground_surface_sem_index;
      new_geometry["boundaries"] = ground_surface;
      co.value()["geometry"].push_back(new_geometry);
    }
  }
}

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
            auto surface = json::array();
            auto surface_sem_index = json::array();
            for (auto &g:co.value()["geometry"]) {
                // get the ground surface from lod0.2 geometry.
                if (g["lod"] == "0.2") {
                    surface = g["boundaries"];
                    surface_sem_index = g["semantics"]["values"];
                    break;
                }
            }
            json first_geo;
            for (auto &g: co.value()["geometry"]) {
                if (g["type"] == "Solid" | g["type"] == "MultiSurface") {
                    first_geo = g;
                    break;
                }
            }
            // record the height of the vertices and their index.
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
            int roof_index = roof_point(height_index);
            // find the surface where the roof point is located.
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
                //// calculate the area of the surface where the roof point is located.
                //// select the surface with the largest area as the main roof surface.
                double max_area = -INFINITY;
                int max_index = -1;
                for (int i = 0; i < surface_index.size(); i++) {
                    int shell_index = surface_index[i].first;
                    int ring_index = surface_index[i].second;
                    Polygon_2 surface_2d;
                    for (auto &vi_list: first_geo["boundaries"][shell_index][ring_index]) {
                        for (auto &vi: vi_list) {
                            std::vector<int> v = j["vertices"][vi.get<int>()];
                            double x = (v[0] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][0].get<double>();
                            double y = (v[1] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][1].get<double>();
                            surface_2d.push_back(K::Point_2(x, y));
                        }
                    }
                    double area = fabs(surface_2d.area());
                    if (area > max_area) {
                        max_area = area;
                        max_index = i;
                    }
                }
                auto main_roof = first_geo["boundaries"][surface_index[max_index].first][surface_index[max_index].second];

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
            }
        }
    }
}