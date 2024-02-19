#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct FaceInfo {
    bool interior, processed;
    FaceInfo() {
        processed = false;
        interior = false;
    }
};
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;

const std::string input_file = "/Users/Acer/Documents/GitHub/Geo1004/Geo1004/hw01/NL.IMBAG.Pand.0503100000032914-0.obj"; // the faculty building
const std::string output_file = "/Users/Acer/Documents/GitHub/Geo1004/Geo1004/hw01/faculty.obj";

struct Vertex {
    double x, y, z;
};

struct Face {
    std::list<int> boundary; // ids of vertices vector
    Kernel::Plane_3 best_plane;
    Triangulation triangulation;
};

int main(int argc, const char * argv[]) {

    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    // Read file
    std::ifstream input_stream;
    input_stream.open(input_file);
    if (input_stream.is_open()) {
        std::string line;

        // Parse line by line
        while (getline(input_stream, line)) {
//      std::cout << line << std::endl;

            std::istringstream line_stream(line);
            char line_type;
            line_stream >> line_type;

            // Vertex
            if (line_type == 'v') {
                vertices.emplace_back();
                double x, y, z;
                line_stream >> x >> y >> z;
                vertices.back().x = x;
                vertices.back().y = y;
                vertices.back().z = z;
            }

            // Face
            if (line_type == 'f') {
                faces.emplace_back();
                int v;
                while (!line_stream.eof()) {
                    line_stream >> v;
                    faces.back().boundary.emplace_back(v-1);
                }
            }

        }
    }

    // Print vertices
//  int i = 0;
//  for (auto const &vertex: vertices) {
//    std::cout << "Vertex " << i++ << ": " << "(" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
//  }

    // Print faces
//  i = 0;
//  for (auto const &face: faces) {
//    std::cout << "Face " << i++ << ": ";
//    for (auto const &vertex: face.boundary) std::cout << " " << vertex;
//    std::cout << std::endl;
//  }

    // Best fitting planes (to do)

    // Triangulate faces (to do)

    // Label triangulation (to do)

    // Export triangles (to do)

    return 0;
}
