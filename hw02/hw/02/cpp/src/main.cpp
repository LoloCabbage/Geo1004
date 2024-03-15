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


int main(int argc, const char * argv[]) {
  //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
  const char* filename = (argc > 1) ? argv[1] : "../../data/tudcampus.city.json";
  std::cout << "Processing: " << filename << std::endl;
  std::ifstream input(filename);
  json j;
  input >> j; //-- store the content of the file in a nlohmann::json object
  input.close();

  generate_lod02(j);

  //-- write to disk the modified city model (out.city.json)
  std::ofstream o("out.city.json");
  o << j.dump(2) << std::endl;
  o.close();

  return 0;
}


// This method get the maximum height and the minimum height of the roof surfaces, and
// use them to calculate the 70 percent roof height by: (max_rh - min_rh) * 0.7 + min_rh.
std::vector<double> roofheight_70p(const json& geometry, const json& j) {
  // initialize max roof height and min roof height.
  double max_roofheight = -INFINITY;
  double min_roofheight = INFINITY;

  // get vertices and transform from the json file.
  auto& vertices = j["vertices"];
  auto& transform = j["transform"];

  if (geometry["type"] == "Solid") {
    // start traverse every surface in the boundaries.
    for (int i = 0; i < geometry["boundaries"].size(); i++) {
      for (int j = 0; j < geometry["boundaries"][i].size(); j++) {

        // get the minimum z value of all vertices from all roof surfaces.
        int sem_index = geometry["semantics"]["values"][i][j];
        if (geometry["semantics"]["surfaces"][sem_index]["type"] == "RoofSurface"){
          auto& roofsurface = geometry["boundaries"][i][j];

          for (auto& rs_part : roofsurface) {
            for (auto& index : rs_part) {
              int vertex_index = index.get<int>();
              auto& vertex = vertices[vertex_index];
              int vertex_z = vertex[2].get<int>();
              double z = (vertex_z * transform["scale"][2].get<double>()) + transform["translate"][2].get<double>();
              if (z > max_roofheight) { max_roofheight = z; }
              if (z < min_roofheight) { min_roofheight = z; }
            }
          }
        }
      }
    }
  }

  if (geometry["type"] == "MultiSurface") {
    for (int i = 0; i < geometry["boundaries"].size(); i++) {
      int sem_index = geometry["semantics"]["values"][i];
      if (geometry["semantics"]["surfaces"][sem_index]["type"] == "RoofSurface") {
        auto& roofsurface = geometry["boundaries"][i];

        for (auto& rs_part : roofsurface) {
          for (auto &index: rs_part) {
            int vertex_index = index.get<int>();
            auto &vertex = vertices[vertex_index];
            int vertex_z = vertex[2].get<int>();
            double z = (vertex_z * transform["scale"][2].get<double>()) + transform["translate"][2].get<double>();
            if (z > max_roofheight) { max_roofheight = z; }
            if (z < min_roofheight) { min_roofheight = z; }
            }
        }
      }
    }
  }

  double roofheight_70p = (max_roofheight - min_roofheight) * 0.7 + min_roofheight;
  return {max_roofheight, min_roofheight, roofheight_70p};
}


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

      // get the 70 percent roof height.
      std::vector<double> roof_height_figures = roofheight_70p(first_geo, j);
      double roof_height_70p = roof_height_figures[2];
//      std::cout << co.key() << " max_rh, min_rh, rh_70p: ";
//      for (double& figure: roof_height_figures) {
//          std::cout << figure << "; ";
//      }
//      std::cout << std::endl;

      // get ground surfaces and their semantic index.
//      if (first_geo["type"] == "Solid") {
//        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
//          for (int k = 0; k < first_geo["boundaries"][i].size(); k++) {
//            int sem_index = first_geo["semantics"]["values"][i][k];
//            if (first_geo["semantics"]["surfaces"][sem_index]["type"] == "GroundSurface") {
//              ground_surface.push_back(first_geo["boundaries"][i][k]);
//              ground_surface_sem_index.push_back(0);
//            }
//          }
//        }
//      } else {
//          for (int i = 0; i < first_geo["boundaries"].size(); i++) {
//            int sem_index = first_geo["semantics"]["values"][i];
//            if (first_geo["semantics"]["surfaces"][sem_index]["type"] == "GroundSurface") {
//              ground_surface.push_back(first_geo["boundaries"][i]);
//              ground_surface_sem_index.push_back(0);
//            }
//          }
//      }

      // get ground surfaces using mean z-value, instead of referring to semantics.
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
            if (surface_mean_z < z_mean && surface_2d_area > area) {
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

      // use ground surface to generate the "roof surfaces" of lod0.2
//      for (auto& gs : new_geometry["boundaries"]) {
//        auto new_gs = json::array();
//        auto& vertices = j["vertices"];
//        auto& transform = j["transform"];
//
//        // use loop here to handle surface with interior boundaries.
//        for (auto& index_list : gs) {
//          auto gs_lift = json::array();
//          for (auto& index: index_list) {
//            int vertex_index = index.get<int>();
//            std::vector<int> vi = vertices[vertex_index];
//            double z_lift = roof_height_70p;
//            int z_int = std::round((z_lift - transform["translate"][2].get<double>()) / transform["scale"][2].get<double>());
//            vertices.push_back(json::array({vi[0], vi[1], z_int}));
//            gs_lift.push_back(vertices.size() - 1);
//          }
//          new_gs.push_back(gs_lift);
//        }
//
//        new_geometry["boundaries"].push_back(new_gs);
//        new_geometry["semantics"]["values"].push_back(1);
//      }

      co.value()["geometry"].push_back(new_geometry);
    }
  }
}