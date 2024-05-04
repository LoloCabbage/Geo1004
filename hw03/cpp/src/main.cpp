#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include "json.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Iso_cuboid_3 Cuboid;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Plane_3<K> Plane;
typedef K::Point_3 Point;
typedef nlohmann::json json;


struct Material {
    std::string name;
    std::vector<double> Kd;
    std::vector<double> Ks;
    int Ns = -999;
    double d = -999;
};

struct Vertex {
    double x, y, z;
};

struct n_Vertex {
    double nx, ny, nz;
};

struct Face {
    std::vector<int> boundary; // ids of vertices vector
    std::vector<int> normalIndices; // ids of normals vector
};

struct Group {
    std::string name;
    std::vector<Face> faces;
    Material material;
};

struct ObjInfo {
    std::vector<Vertex> vertices;
    std::vector<n_Vertex> normals;
    std::vector<Group> allGroups;
};

struct Voxel {
    unsigned int row_x_num, row_y_num, row_z_num, value, group;
    double voxel_size;
    double min_x, min_y, min_z, max_x, max_y, max_z;

    Voxel(unsigned int x, unsigned int y, unsigned int z,
          double grid_minx, double grid_miny, double grid_minz,
          unsigned int default_value, double resolution, unsigned int group_label)
    {
      value = default_value;
      group = group_label;
      row_x_num = x;
      row_y_num = y;
      row_z_num = z;
      voxel_size = resolution;
      min_x = x * voxel_size + grid_minx;
      min_y = y * voxel_size + grid_miny;
      min_z = z * voxel_size + grid_minz;
      max_x = min_x + voxel_size;
      max_y = min_y + voxel_size;
      max_z = min_z + voxel_size;
    }
};

struct VoxelSurface {
    std::vector<Point> vertices; // Stores unique vertices
    std::vector<std::vector<size_t>> faces; // Indices to vertices for each face
    std::vector<Plane> planes; // CGAL planes for each face
//    std::vector<K::Vector_3> normals; // Normal for each face
};

struct Room {
    std::map<unsigned int, VoxelSurface> roomSurfaces; // Maps room value to its geometry
    VoxelSurface outerEnvelope; // Separate geometry data for the outer envelope
};

ObjInfo read_obj(const std::string& obj, const std::string& mtl) {

  std::vector<Vertex> vertices;
  std::vector<n_Vertex> normals;
  std::vector<Group> allGroups;

  // Read obj file
  std::ifstream input_stream_obj;
  std::ifstream input_stream_mtl;
  input_stream_obj.open(obj);
  input_stream_mtl.open(mtl);
  if (input_stream_obj.is_open() && input_stream_mtl.is_open()) {
    std::string line;
    std::string material;
    Group currentGroup;
    Material current_mtl;
    bool isGroupActive = false; // flag to check if a group is active
    bool isMtlActive = false;

    // get all material info
    std::vector<Material> all_materials;
    while (getline(input_stream_mtl, line)) {
      std::istringstream line_stream(line);
      std::string line_type;
      line_stream >> line_type;

      if (isMtlActive) {
        if (line_type == "Kd") {
          double Kd1, Kd2, Kd3;
          line_stream >> Kd1 >> Kd2 >> Kd3;
          current_mtl.Kd.emplace_back(Kd1);
          current_mtl.Kd.emplace_back(Kd2);
          current_mtl.Kd.emplace_back(Kd3);
        }

        if (line_type == "Ks") {
          double Ks1, Ks2, Ks3;
          line_stream >> Ks1 >> Ks2 >> Ks3;
          current_mtl.Ks.emplace_back(Ks1);
          current_mtl.Ks.emplace_back(Ks2);
          current_mtl.Ks.emplace_back(Ks3);
        }

        if (line_type == "Ns") {
          line_stream >> current_mtl.Ns;
        }

        if (line_type == "d") {
          line_stream >> current_mtl.d;
        }
      }

      if (line_type == "newmtl") {
        if (isMtlActive) {
          all_materials.emplace_back(current_mtl);
          current_mtl.Kd.clear();
          current_mtl.Ks.clear();
        }
        line_stream >> current_mtl.name;
        isMtlActive = true;
      }
    }
    if (isMtlActive) {
      all_materials.emplace_back(current_mtl);
    }

    // Parse line by line
    while (getline(input_stream_obj, line)) {

      std::istringstream line_stream(line);
      std::string line_type;
      line_stream >> line_type;

      // Vertex
      if (line_type == "v") {
        vertices.emplace_back();
        double x, y, z;
        line_stream >> x >> y >> z;
        vertices.back().x = x;
        vertices.back().y = y;
        vertices.back().z = z;
//        std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
      }

      // Normal
      if (line_type == "vn") {
        normals.emplace_back();
        double nx, ny, nz;
        line_stream >> nx >> ny >> nz;
        normals.back().nx = nx;
        normals.back().ny = ny;
        normals.back().nz = nz;
//        std::cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;
      }

      // Face
      if (line_type == "f") {
        Face face;
        std::string vertex_info;
        while (line_stream >> vertex_info) {
          std::replace(vertex_info.begin(), vertex_info.end(), '/', ' ');//replace // with /
          std::istringstream vertex_stream(vertex_info);
          int v, vn;
          vertex_stream >> v; // read v value
          vertex_stream.ignore(std::numeric_limits<std::streamsize>::max(), ' '); // ignore the next character
          vertex_stream.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
          vertex_stream >> vn; // read vn value
          face.boundary.push_back(v);
          face.normalIndices.push_back(vn);
        }
        if (isGroupActive) { // if a group is active, add the face to the group
          currentGroup.faces.push_back(face);
        }
      }

      // Material
      if (line_type == "usemtl") {
        Material current_m;
        std::string mtl_name;
        line_stream >> mtl_name;
        for (Material& m : all_materials) {
          if (mtl_name == m.name) { current_m = m; }
        }
        if (isGroupActive) { currentGroup.material = current_m; }
      }

      // start a new Group
      if (line_type == "g") {
        if (isGroupActive) {
          allGroups.push_back(currentGroup);
          currentGroup.faces.clear(); // prepare for the next group
        }
        line_stream >> currentGroup.name;
        isGroupActive = true;
      }
    }
    if (isGroupActive) {
      allGroups.push_back(currentGroup);
    }
    // Print out the number of faces in each group
//    for (const auto& group : allGroups) {
//      std::cout << "Group " << group.name << " has " << group.faces.size() << " faces." << std::endl;
//    }
  }
  else{
    throw std::runtime_error("Unable to open obj file");
  }

  ObjInfo obj_info;
  obj_info.vertices = vertices;
  obj_info.normals = normals;
  obj_info.allGroups = allGroups;
  return obj_info;
}

struct VoxelGrid {
  std::vector<Voxel> voxels; // stores the voxel grid
  unsigned int max_x, max_y, max_z; // max number of voxels in x,y,z directions
  double grid_res;
  double coord_minx, coord_miny, coord_minz; // min/max xyz coordinates
  std::map<int, unsigned int> roomVoxelCounts;

  VoxelGrid(double bbox_minx, double bbox_miny, double bbox_minz,
            double bbox_maxx, double bbox_maxy, double bbox_maxz, double voxel_size)
  {
    // calculate the number of voxels in x,y,z directions
    // make sure that the voxel grid covers the entire bounding box
    unsigned int x = ceil((bbox_maxx - bbox_minx) / voxel_size);
    unsigned int y = ceil((bbox_maxy - bbox_miny) / voxel_size);
    unsigned int z = ceil((bbox_maxz - bbox_minz) / voxel_size);
    max_x = x+10;
    max_y = y+10;
    max_z = z+10;
    coord_minx = bbox_minx-5*voxel_size;
    coord_miny = bbox_miny-5*voxel_size;
    coord_minz = bbox_minz-5*voxel_size;
    grid_res = voxel_size;
    for (unsigned int i = 0; i < max_z; ++i) {
      for (unsigned int j = 0; j < max_y; ++j) {
        for (unsigned int k = 0; k < max_x; ++k) {
          Voxel v(k, j, i,
                  coord_minx, coord_miny, coord_minz,
                  0, voxel_size, 0);
          voxels.push_back(v);
        }
      }
    }
  }

  // x,y,z are the row numbers of the voxel
  Voxel &operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z)
  { // allows to read and write the voxel value
    assert(x >= 0 && x < max_x); // check if x,y,z is within the range
    assert(y >= 0 && y < max_y);
    assert(z >= 0 && z < max_z);
    return voxels[x + y*max_x + z*max_x*max_y];
  }

  Voxel operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z)
  const { // read only version of the above function
    assert(x >= 0 && x < max_x);
    assert(y >= 0 && y < max_y);
    assert(z >= 0 && z < max_z);
    return voxels[x + y*max_x + z*max_x*max_y];
  }

  void triangle_voxel_intersection(const ObjInfo& obj_info, VoxelGrid& voxel_grid) {
    int group_label = 0;
    for (const auto &group: obj_info.allGroups) {
      group_label++;
      for (const auto &face: group.faces) {
        if (face.boundary.size() == 3) {
          Kernel::Point_3 p1(obj_info.vertices[face.boundary[0] - 1].x,
                             obj_info.vertices[face.boundary[0] - 1].y,
                             obj_info.vertices[face.boundary[0] - 1].z);
          Kernel::Point_3 p2(obj_info.vertices[face.boundary[1] - 1].x,
                             obj_info.vertices[face.boundary[1] - 1].y,
                             obj_info.vertices[face.boundary[1] - 1].z);
          Kernel::Point_3 p3(obj_info.vertices[face.boundary[2] - 1].x,
                             obj_info.vertices[face.boundary[2] - 1].y,
                             obj_info.vertices[face.boundary[2] - 1].z);
          Kernel::Triangle_3 triangle(p1, p2, p3);
          std::vector<CGAL::Point_3<CGAL::Epick>> tri_points = {p1, p2, p3};
          // get the bounding box of the triangle
          Kernel::Iso_cuboid_3 t_bbox = CGAL::bounding_box(tri_points.begin(), tri_points.end());
          Kernel::Point_3 t_bottom_left = t_bbox.min();
          Kernel::Point_3 t_top_right = t_bbox.max();
          // get the voxel indices of the bounding box
          unsigned int min_x = floor((t_bottom_left.x() - voxel_grid.coord_minx) / voxel_grid.grid_res);
          unsigned int min_y = floor((t_bottom_left.y() - voxel_grid.coord_miny) / voxel_grid.grid_res);
          unsigned int min_z = floor((t_bottom_left.z() - voxel_grid.coord_minz) / voxel_grid.grid_res);
          unsigned int max_x = ceil((t_top_right.x() - voxel_grid.coord_minx) / voxel_grid.grid_res);
          unsigned int max_y = ceil((t_top_right.y() - voxel_grid.coord_miny) / voxel_grid.grid_res);
          unsigned int max_z = ceil((t_top_right.z() - voxel_grid.coord_minz) / voxel_grid.grid_res);

          // iterate over the voxels in the expanded bounding box
          for (unsigned int i = min_z-1; i < max_z+1; ++i) {
            for (unsigned int j = min_y-1; j < max_y+1; ++j) {
              for (unsigned int k = min_x-1; k < max_x+1; ++k) {
                Voxel& voxel = voxel_grid(k, j, i);
                if (voxel.value == 1) {
                  continue; // if the voxel is already intersecting with a face, skip it
                }
                auto cube = Cuboid(voxel.min_x, voxel.min_y, voxel.min_z,
                                   voxel.max_x, voxel.max_y,voxel.max_z);
                //std::cout << cube << std::endl;
                if (CGAL::do_intersect(triangle, cube)) {
                  //std::cout << "Face is intersecting with a voxel grid" << k << " " << j << " " << i << std::endl;
                  voxel.value = 1;
                  voxel.group = group_label;
                }
              }
            }
          }
        }
        else{
          std::cout << "Face is not a triangle" << std::endl;
        }
      }
    }
    //std::cout << "Final group label: " << group_label << std::endl;
  }

  bool isBoundaryVoxel(const unsigned int &x, const unsigned int &y, const unsigned int &z)
  const {
    return x == 0 || y == 0 || z == 0 || x == max_x - 1 || y == max_y - 1 || z == max_z - 1;
  }// check if the voxel is on the boundary of the voxel grid

  bool isValidVoxel(const unsigned int &x, const unsigned int &y, const unsigned int &z)
  const {
    return x >= 0 && x < max_x && y >= 0 && y < max_y && z >= 0 && z < max_z;
  } // check if the voxel is within the voxel grid

  void fillExterior(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
    std::list<std::tuple<unsigned int, unsigned int, unsigned int>> queue;
    queue.emplace_back(x, y, z); // add the current voxel to the queue

    while(!queue.empty()) {
      auto voxel_xyz = queue.front(); // get the first element in the queue
      queue.pop_front(); // remove the first element from the queue
      unsigned int current_x = std::get<0>(voxel_xyz); // get the x,y,z values of the voxel
      unsigned int current_y = std::get<1>(voxel_xyz);
      unsigned int current_z = std::get<2>(voxel_xyz);
      unsigned int index = current_x + current_y*max_x + current_z*max_x*max_y;
      if (isValidVoxel(current_x, current_y, current_z) && voxels[index].value == 0){
        voxels[index].value = 2; // mark the voxel as exterior
        queue.emplace_back(current_x + 1, current_y, current_z); // add the neighbors to the queue
        queue.emplace_back(current_x - 1, current_y, current_z);
        queue.emplace_back(current_x, current_y + 1, current_z);
        queue.emplace_back(current_x, current_y - 1, current_z);
        queue.emplace_back(current_x, current_y, current_z + 1);
        queue.emplace_back(current_x, current_y, current_z - 1);
      }
    }
  }

  void markExterior() {
    for (unsigned int i = 0; i < max_z; ++i) {
      for (unsigned int j = 0; j < max_y; ++j) {
        for (unsigned int k = 0; k < max_x; ++k) {
          if (isBoundaryVoxel(k, j, i) && voxels[k + j*max_x + i*max_x*max_y].value == 0) {
            fillExterior(k, j, i);
          }
        }
      }
    }
  }

  void fillRooms(const unsigned int &x, const unsigned int &y, const unsigned int &z, int &roomID) {
    std::list<std::tuple<unsigned int, unsigned int, unsigned int>> queue;
    queue.emplace_back(x, y, z);

    while(!queue.empty()) {
      auto voxel_xyz = queue.front();
      queue.pop_front();
      unsigned int current_x = std::get<0>(voxel_xyz);
      unsigned int current_y = std::get<1>(voxel_xyz);
      unsigned int current_z = std::get<2>(voxel_xyz);
      unsigned int index = current_x + current_y*max_x + current_z*max_x*max_y;
      if (isValidVoxel(current_x, current_y, current_z) && voxels[index].value == 0){
        voxels[index].value = roomID;
        roomVoxelCounts[roomID]++;
        queue.emplace_back(current_x + 1, current_y, current_z);
        queue.emplace_back(current_x - 1, current_y, current_z);
        queue.emplace_back(current_x, current_y + 1, current_z);
        queue.emplace_back(current_x, current_y - 1, current_z);
        queue.emplace_back(current_x, current_y, current_z + 1);
        queue.emplace_back(current_x, current_y, current_z - 1);
      }
    }
  }

  void markRooms() {
    int roomID = 3;
    for (unsigned int i = 0; i < max_z; ++i) {
      for (unsigned int j = 0; j < max_y; ++j) {
        for (unsigned int k = 0; k < max_x; ++k) {
          if (voxels[k + j*max_x + i*max_x*max_y].value == 0) {
            fillRooms(k, j, i, roomID);
            if (checkValidRoom(roomID)) roomID++;
            else {
              roomVoxelCounts.erase(roomID);
              for (unsigned int o = 0; o < max_z; ++o) {
                for (unsigned int p = 0; p < max_y; ++p) {
                  for (unsigned int q = 0; q < max_x; ++q) {
                    unsigned int index = q + p*max_x + o*max_x*max_y;
                    Voxel& voxel = voxels[index];
                    if (voxel.value == roomID) voxel.value = 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  static double distance(const Point& p1, const Point& p2) {
    double dx = p1.x() - p2.x();
    double dy = p1.y() - p2.y();
    double dz = p1.z() - p2.z();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  bool checkValidRoom(const int& roomID) {
    std::vector<Point> RoomPoints;
    std::array<Point, 8> ObbPoints;
    for (unsigned int i = 0; i < max_z; ++i) {
      for (unsigned int j = 0; j < max_y; ++j) {
        for (unsigned int k = 0; k < max_x; ++k) {
          unsigned int index = k + j*max_x + i*max_x*max_y;
          Voxel& voxel = voxels[index];
          if (voxel.value == roomID) {
            Point p0(voxel.min_x, voxel.min_y, voxel.min_z);
            Point p1(voxel.max_x, voxel.min_y, voxel.min_z);
            Point p2(voxel.max_x, voxel.max_y, voxel.min_z);
            Point p3(voxel.min_x, voxel.max_y, voxel.min_z);
            Point p4(voxel.min_x, voxel.max_y, voxel.max_z);
            Point p5(voxel.min_x, voxel.min_y, voxel.max_z);
            Point p6(voxel.max_x, voxel.min_y, voxel.max_z);
            Point p7(voxel.max_x, voxel.max_y, voxel.max_z);
            RoomPoints.push_back(p0);
            RoomPoints.push_back(p1);
            RoomPoints.push_back(p2);
            RoomPoints.push_back(p3);
            RoomPoints.push_back(p4);
            RoomPoints.push_back(p5);
            RoomPoints.push_back(p6);
            RoomPoints.push_back(p7);
          }
        }
      }
    }
    CGAL::oriented_bounding_box(RoomPoints, ObbPoints);
    double bbox_x = distance(ObbPoints[1], ObbPoints[0]);
    double bbox_y = distance(ObbPoints[2], ObbPoints[1]);
    double bbox_z = distance(ObbPoints[4], ObbPoints[3]);
    if (bbox_x < 1 || bbox_y < 1 || bbox_z < 1) return false;
    else return true;
  }
};

VoxelGrid create_voxel_grid(const ObjInfo& obj_info, double voxel_size) {
  // read the vertices and get the bounding box
  std::vector<Kernel::Point_3> points;
  for (const auto& vertex : obj_info.vertices) {
    points.emplace_back(vertex.x, vertex.y, vertex.z);
  }
  Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
  Kernel::Point_3 bottom_left = bbox.min();
  Kernel::Point_3 top_right = bbox.max();
  double min_x = bottom_left.x();
  double min_y = bottom_left.y();
  double min_z = bottom_left.z();
  double max_x = top_right.x();
  double max_y = top_right.y();
  double max_z = top_right.z();

  VoxelGrid voxel_grid(min_x, min_y, min_z,
                       max_x, max_y, max_z, voxel_size);

  return voxel_grid;
}

// Write voxels into obj file for visualization.
void WriteVoxelsToOBJ(const VoxelGrid& voxelGrid, const std::string& filename, const std::string& voxeltype) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return;
  }

  // OBJ files are 1-indexed
  unsigned int vertexIndex = 1;
  for (const auto& voxel : voxelGrid.voxels) {
    bool condition = false;
    if (voxeltype == "buildingpart" && voxel.value == 1) {condition = true;}
    else if (voxeltype == "exterior" && voxel.value == 2) {condition = true;}
    else if (voxeltype == "rooms" && voxel.value != 0 && voxel.value != 1 && voxel.value != 2) {condition = true;}
    if (condition) {
      // Define vertices for the current voxel
      std::vector<std::tuple<double, double, double>> vertices = {
              {voxel.min_x, voxel.min_y, voxel.min_z},
              {voxel.max_x, voxel.min_y, voxel.min_z},
              {voxel.max_x, voxel.max_y, voxel.min_z},
              {voxel.min_x, voxel.max_y, voxel.min_z},
              {voxel.min_x, voxel.min_y, voxel.max_z},
              {voxel.max_x, voxel.min_y, voxel.max_z},
              {voxel.max_x, voxel.max_y, voxel.max_z},
              {voxel.min_x, voxel.max_y, voxel.max_z},
      };

      // Write vertices to file
      for (const auto& vertex : vertices) {
        file << "v " << std::get<0>(vertex) << " " << std::get<1>(vertex) << " " << std::get<2>(vertex) << std::endl;
      }

      // Define faces using vertex indices (note: OBJ format uses 1-based indexing)
      std::vector<std::vector<unsigned int>> faces = {
        {vertexIndex, vertexIndex + 1, vertexIndex + 2, vertexIndex + 3},
        {vertexIndex + 4, vertexIndex + 5, vertexIndex + 6, vertexIndex + 7},
        {vertexIndex, vertexIndex + 1, vertexIndex + 5, vertexIndex + 4},
        {vertexIndex + 2, vertexIndex + 3, vertexIndex + 7, vertexIndex + 6},
        {vertexIndex + 1, vertexIndex + 2, vertexIndex + 6, vertexIndex + 5},
        {vertexIndex, vertexIndex + 3, vertexIndex + 7, vertexIndex + 4}
      };

      // Write faces to file
      for (const auto& face : faces) {
        file << "f";
        for (unsigned int index : face) {
          file << " " << index;
        }
        file << std::endl;
      }

      vertexIndex += 8; // Move index for the next voxel
    }
  }

  file.close();
}

// Adds a vertex to Room if it doesn't already exist and returns its index
size_t addVertex(VoxelSurface& geomData, const Point& point) {
  auto it = std::find(geomData.vertices.begin(), geomData.vertices.end(), point);
  if (it != geomData.vertices.end()) {
    return std::distance(geomData.vertices.begin(), it); // Found existing point
  } else {
    geomData.vertices.push_back(point); // Add new point
    return geomData.vertices.size() - 1; // Index of new point
  }
}

//add the identified surface to the corresponding room or outer surface
void add_surface(VoxelSurface& geomData, const VoxelGrid& grid, unsigned int k, unsigned int j, unsigned int i, int dx, int dy, int dz) {
  // Calculate the base position of the voxel
  double baseX = grid.coord_minx + k * grid.grid_res;
  double baseY = grid.coord_miny + j * grid.grid_res;
  double baseZ = grid.coord_minz + i * grid.grid_res;

  // Vertices of the face
  std::vector<Point> Points(4);

  // Create a plane from three non-collinear points of the face
  Plane plane = Plane(Points[0], Points[1], Points[2]);
//    std::vector<Vector> normals;

  if (dx != 0) { //Left or Right face
    double x = baseX + (dx > 0 ? grid.grid_res : 0);
    Points[0] = Point(x, baseY, baseZ);
    Points[1] = Point(x, baseY + grid.grid_res, baseZ);
    Points[2] = Point(x, baseY + grid.grid_res, baseZ + grid.grid_res);
    Points[3] = Point(x, baseY, baseZ + grid.grid_res);
    if (dx < 0) std::reverse(Points.begin(), Points.end()); // Reverse for inward normal
  }
  else if (dy != 0) { //bottom or top face
    double y = baseY + (dy > 0 ? grid.grid_res : 0);
    Points[0] = Point(baseX, y, baseZ);
    Points[1] = Point(baseX + grid.grid_res, y, baseZ);
    Points[2] = Point(baseX + grid.grid_res, y, baseZ + grid.grid_res);
    Points[3] = Point(baseX, y, baseZ + grid.grid_res);
    if (dy > 0) std::reverse(Points.begin(), Points.end());
  }

  else if (dz != 0) { // Front or Back face
    double z = baseZ + (dz > 0 ? grid.grid_res : 0);
    Points[0] = Point(baseX, baseY, z);
    Points[1] = Point(baseX, baseY + grid.grid_res, z);
    Points[2] = Point(baseX + grid.grid_res, baseY + grid.grid_res, z);
    Points[3] = Point(baseX + grid.grid_res, baseY, z);
    if (dz > 0) std::reverse(Points.begin(), Points.end());
  }

//  // Manually determine the outward normal vector based on voxel neighborhood
//  Vector normalVector(dx * grid.grid_res, dy * grid.grid_res, dz * grid.grid_res);
//  geomData.normals.push_back(normalVector);

  // Adding vertices to VoxelSurface and constructing the face
  std::vector<size_t> faceIndices;
  for (const auto& point : Points) {
      size_t idx = addVertex(geomData, point);
      faceIndices.push_back(idx);
  }
  geomData.faces.push_back(faceIndices);
  geomData.planes.push_back(plane);
//    geomData.normals.push_back(normalVector);
}

//Identifies the surface between each building part's voxel and its neighbor with different values, in 6 directions
void surface_extraction(const VoxelGrid& grid, Room& roomGeom) {
  for (unsigned int i = 0; i < grid.max_z; ++i) {
    for (unsigned int j = 0; j < grid.max_y; ++j) {
      for (unsigned int k = 0; k < grid.max_x; ++k) {
        unsigned int currentIndex = k + j * grid.max_x + i * grid.max_x * grid.max_y;
        unsigned int currentValue = grid.voxels[currentIndex].value;

        if (currentValue == 1) { // Building part voxel
          std::vector<std::tuple<int, int, int>> neighborDirections = {
                  {1, 0, 0}, {-1, 0, 0},
                  {0, 1, 0}, {0, -1, 0},
                  {0, 0, 1}, {0, 0, -1}};
          for (const auto& dir : neighborDirections) {
            int dx = std::get<0>(dir), dy = std::get<1>(dir), dz = std::get<2>(dir);
            unsigned int nx = k + dx, ny = j + dy, nz = i + dz;

            if (nx < grid.max_x && ny < grid.max_y && nz < grid.max_z) {
              unsigned int neighborIndex = nx + ny * grid.max_x + nz * grid.max_x * grid.max_y;
              unsigned int neighborValue = grid.voxels[neighborIndex].value;

              if (neighborValue == 2) {
                // Outer envelope surface between building part and exterior
                add_surface(roomGeom.outerEnvelope, grid, k, j, i, dx, dy, dz);
              } else if (neighborValue >= 3) {
                // If the voxel is at the boundary of different rooms, this will handle creating surfaces for each
                add_surface(roomGeom.roomSurfaces[neighborValue], grid, k, j, i, dx, dy, dz);
              }
            }
          }
        }
      }
    }
  }
}

// write a single room in the struct to a single obj file
void writeGeometryToOBJ(const VoxelSurface& geomData, const std::string& filename) {
  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // Write vertices
  for (const auto& vertex : geomData.vertices) {
    outFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
  }

  // Write faces
  for (const auto& face : geomData.faces) {
    outFile << "f";
    for (auto index : face) {
        outFile << " " << (index + 1); // OBJ format uses 1-based indexing
    }
    outFile << "\n";
  }

  outFile.close();
}

void exportRoomGeometries(const Room& roomGeom, const std::string& basePath) {
  // Write the outer envelope geometry
  std::string envelopeFilename = basePath + "outer_envelope.obj";
  writeGeometryToOBJ(roomGeom.outerEnvelope, envelopeFilename);

  // Write each room's geometry
  for (const auto& room : roomGeom.roomSurfaces) {
    unsigned int roomID = room.first;
    const VoxelSurface& geomData = room.second;

    std::string roomFilename = basePath + "room_" + std::to_string(roomID) + ".obj";
    writeGeometryToOBJ(geomData, roomFilename);
  }
}

void writeCityJson(const Room& roomGeom, const VoxelGrid& voxel_grid) {
  json j;
  double minx = voxel_grid.coord_minx;
  double miny = voxel_grid.coord_miny;
  double minz = voxel_grid.coord_minz;
  double scale = voxel_grid.grid_res;
  j["type"] = "CityJSON";
  j["version"] = "2.0";
  j["transform"] = nlohmann::json::object();
  j["transform"]["scale"] = nlohmann::json::array({scale, scale, scale});
  j["transform"]["translate"] = nlohmann::json::array({minx, miny, minz});
  j["CityObjects"] = nlohmann::json::object();
  j["vertices"] = json::array();
  j["CityObjects"]["BuildingParent"] = {{"attributes", json::object()},
                                        {"type", "Building"},
                                        {"geometry", json::array()}};
  j["CityObjects"]["BuildingParent"]["children"] = {"OuterEnvelope"};

  j["CityObjects"]["OuterEnvelope"] = {{"attributes", json::object()}, {"geometry", json::array()},
                                       {"type", "BuildingPart"}, {"parents", {"BuildingParent"}}};
  j["CityObjects"]["OuterEnvelope"]["geometry"].push_back({{"lod", "2"},
                                                           {"type", "MultiSurface"},
                                                           {"boundaries", json::array()}});
  for (auto& face : roomGeom.outerEnvelope.faces) {
    unsigned int current_vertices = j["vertices"].size();
    std::vector<size_t> new_face;
    for (auto& num : face) {
      new_face.push_back(num + current_vertices);
    }
    j["CityObjects"]["OuterEnvelope"]["geometry"][0]["boundaries"].push_back({new_face});
  }
  for (auto& vertex : roomGeom.outerEnvelope.vertices) {
    size_t v0 = round((vertex[0] - minx) / scale);
    size_t v1 = round((vertex[1] - miny) / scale);
    size_t v2 = round((vertex[2] - minz) / scale);
    j["vertices"].push_back({v0, v1, v2});
  }

  for (int i = 0; i < roomGeom.roomSurfaces.size(); i++) {
    std::string roomName = "room" + std::to_string(i);
    j["CityObjects"]["BuildingParent"]["children"].emplace_back(roomName);
    j["CityObjects"][roomName] = {{"attributes", json::object()}, {"geometry", json::array()},
                                  {"type", "BuildingRoom"}, {"parents", {"BuildingParent"}}};
    j["CityObjects"][roomName]["geometry"].push_back({{"lod", "2"},
                                                      {"type", "MultiSurface"},
                                                      {"boundaries", json::array()}});
    for (auto& face : roomGeom.roomSurfaces.at(i+3).faces) {
      unsigned int current_vertices = j["vertices"].size();
      std::vector<size_t> new_face;
      for (auto& num : face) {
        new_face.push_back(num + current_vertices);
      }
      j["CityObjects"][roomName]["geometry"][0]["boundaries"].push_back({new_face});
    }
    for (auto& vertex : roomGeom.roomSurfaces.at(i+3).vertices) {
      size_t v0 = round((vertex[0] - minx) / scale);
      size_t v1 = round((vertex[1] - miny) / scale);
      size_t v2 = round((vertex[2] - minz) / scale);
      j["vertices"].push_back({v0, v1, v2});
    }
  }

  std::string json_string = j.dump(2);
  std::ofstream out_stream("output.json");
  out_stream << json_string;
  out_stream.close();
}


int main(int argc, const char * argv[]) {
  const char* input_obj = (argc > 2) ? argv[1] : "../../data/Wellness_center_Sama.obj";
  const char* input_mtl = (argc > 2) ? argv[2] : "../../data/Wellness_center_Sama.mtl";
  if (argc <= 2) {
    std::cout << "No specified input OBJ or MTL file, default files are applied." << std::endl;
    printf("OBJ Path: %s\n", input_obj);
    printf("MTL Path: %s\n", input_mtl);
  }
  ObjInfo obj_info = read_obj(input_obj, input_mtl);
  std::cout << "Number of groups in obj file: " << obj_info.allGroups.size() << std::endl;
  VoxelGrid voxel_grid = create_voxel_grid(obj_info, 0.2);
  std::cout << "Voxel grid xyz numbers: " << voxel_grid.max_x <<","<<voxel_grid.max_y<<","<<voxel_grid.max_z<<std::endl;
  std::cout << "Start marking voxels." << std::endl;
  voxel_grid.triangle_voxel_intersection(obj_info, voxel_grid);
  voxel_grid.markExterior();
  voxel_grid.markRooms();
  std::cout << "There are " << voxel_grid.roomVoxelCounts.size() << " rooms." << std::endl;
  //WriteVoxelsToOBJ(voxel_grid, "output_voxel.obj", "rooms");
  std::cout << "Voxel marking is finished." << std::endl;

  std::cout << "Start surface extraction." << std::endl;
  // Initialize Room to store surface extractions
  Room Room;

  // Perform surface extraction
  surface_extraction(voxel_grid, Room);
  writeCityJson(Room, voxel_grid);

//  // Define a base path for output OBJ files
//  std::string outputPath = "../../output/";
//
//  // Export extracted geometries to separate OBJ files
//  exportRoomGeometries(Room, outputPath);

  return 0;
}

