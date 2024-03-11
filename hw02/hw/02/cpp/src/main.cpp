/*
+------------------------------------------------------------------------------+
|                                                                              |
|                                 Hugo Ledoux                                  |
|                             h.ledoux@tudelft.nl                              |
|                                  2024-02-21                                  |
|                                                                              |
+------------------------------------------------------------------------------+
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder

using json = nlohmann::json;


int   get_no_roof_surfaces(json &j);
void  list_all_vertices(json& j);
void  visit_roofsurfaces(json &j);
void  generate_lod02(json &j);
void generate_lod12(json &j);

int main(int argc, const char * argv[]) {
  //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
  const char* filename = (argc > 1) ? argv[1] : "../../data/tudcampus.city.json";
  std::cout << "Processing: " << filename << std::endl;
  std::ifstream input(filename);
  json j;
  input >> j; //-- store the content of the file in a nlohmann::json object
  input.close();

  //-- get total number of RoofSurface in the file
  int noroofsurfaces = get_no_roof_surfaces(j);
  std::cout << "Total RoofSurface: " << noroofsurfaces << std::endl;

  // list_all_vertices(j);

  // visit_roofsurfaces(j);

  generate_lod02(j);
  generate_lod12(j);
  //-- print out the number of Buildings in the file
  int nobuildings = 0;
  for (auto& co : j["CityObjects"]) {
    if (co["type"] == "Building") {
      nobuildings += 1;
    }
  }
  std::cout << "There are " << nobuildings << " Buildings in the file" << std::endl;

  //-- print out the number of vertices in the file
  std::cout << "Number of vertices " << j["vertices"].size() << std::endl;

  //-- write to disk the modified city model (out.city.json)
  std::ofstream o("out.city.json");
  o << j.dump(2) << std::endl;
  o.close();

  return 0;
}


// Visit every 'RoofSurface' in the CityJSON model and output its geometry (the arrays of indices)
// Useful to learn to visit the geometry boundaries and at the same time check their semantics.
void visit_roofsurfaces(json &j) {
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] == "Solid") {
        for (int i = 0; i < g["boundaries"].size(); i++) {
          for (int j = 0; j < g["boundaries"][i].size(); j++) {
            int sem_index = g["semantics"]["values"][i][j];
            if (g["semantics"]["surfaces"][sem_index]["type"].get<std::string>().compare("RoofSurface") == 0) {
              std::cout << "RoofSurface: " << g["boundaries"][i][j] << std::endl;
            }
          }
        }
      }
    }
  }
}


// Returns the number of 'RooSurface' in the CityJSON model
int get_no_roof_surfaces(json &j) {
  int total = 0;
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] == "Solid") {
        for (auto& shell : g["semantics"]["values"]) {
          for (auto& s : shell) {
            if (g["semantics"]["surfaces"][s.get<int>()]["type"].get<std::string>().compare("RoofSurface") == 0) {
              total += 1;
            }
          }
        }
      }
    }
  }
  return total;
}


// CityJSON files have their vertices compressed: https://www.cityjson.org/specs/1.1.1/#transform-object
// this function visits all the surfaces and print the (x,y,z) coordinates of each vertex encountered
void list_all_vertices(json& j) {
  for (auto& co : j["CityObjects"].items()) {
    std::cout << "= CityObject: " << co.key() << std::endl;
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] == "Solid") {
        for (auto& shell : g["boundaries"]) {
          for (auto& surface : shell) {
            for (auto& ring : surface) {
              std::cout << "---" << std::endl;
              for (auto& v : ring) { 
                std::vector<int> vi = j["vertices"][v.get<int>()];
                double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
                double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
                double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                std::cout << std::setprecision(2) << std::fixed << v << " (" << x << ", " << y << ", " << z << ")" << std::endl;                
              }
            }
          }
        }
      }
    }
  }
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

      //// get the first geometry which has the type "Solid" or "MultiSurface", use this to generate lod0.2.
      json first_geo;
      for (auto& g : co.value()["geometry"]) {
        if (g["type"] == "Solid" | g["type"] == "MultiSurface") {
          first_geo = g;
          break;
        }
      }

      //// get the 70 percent roof height.
      std::vector<double> roof_height_figures = roofheight_70p(first_geo, j);
      double roof_height_70p = roof_height_figures[2];
      std::cout << co.key() << " max_rh, min_rh, rh_70p: ";
      for (double& figure: roof_height_figures) {
          std::cout << figure << "; ";
      }
      std::cout << std::endl;

      //// get ground surfaces and their semantic index.
      if (first_geo["type"] == "Solid") {
        for (int i = 0; i < first_geo["boundaries"].size(); i++) {
          for (int k = 0; k < first_geo["boundaries"][i].size(); k++) {
            int sem_index = first_geo["semantics"]["values"][i][k];
            if (first_geo["semantics"]["surfaces"][sem_index]["type"] == "GroundSurface") {
              ground_surface.push_back(first_geo["boundaries"][i][k]);
              ground_surface_sem_index.push_back(0);
            }
          }
        }
      } else {
          for (int i = 0; i < first_geo["boundaries"].size(); i++) {
            int sem_index = first_geo["semantics"]["values"][i];
            if (first_geo["semantics"]["surfaces"][sem_index]["type"] == "GroundSurface") {
              ground_surface.push_back(first_geo["boundaries"][i]);
              ground_surface_sem_index.push_back(0);
            }
          }
      }


      //// create new lod0.2 geometry
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
        for (auto& g : co.value()["geometry"]) {
        if (g["lod"] == "1.2") {
            lod12_exist = true;
            break;
        }
        }
        if (lod12_exist) {
        std::cout << "This CityObject already contains LoD1.2 geometry, skip." << std::endl;
        continue;
        }

        if (!co.value()["geometry"].empty()) {
        // initialize ground surface array and ground surface semantic index array.
        auto ground_surface = json::array();
        auto ground_surface_sem_index = json::array();

        //// get the first geometry which has the type "Solid" or "MultiSurface", use this to generate lod1.2.
        json first_geo;
        for (auto& g : co.value()["geometry"]) {
            if (g["type"] == "Solid" | g["type"] == "MultiSurface") {
            first_geo = g;
            break;
            }
        }
