#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

namespace fs = std::filesystem;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; // import ... as kernel
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


const std::string input_file = "/mnt/c/Users/Acer/Documents/GitHub/Geo1004/Geo1004/hw01/NL.IMBAG.Pand.0503100000020110-0.obj"; // the faculty building
const std::string output_file = "/mnt/c/Users/Acer/Documents/GitHub/Geo1004/Geo1004/hw01/faculty.obj";

struct Vertex {
    double x, y, z;
};

struct Face {
    std::list<int> boundary; // ids of vertices vector, boundary = []
    Kernel::Plane_3 best_plane;
    Triangulation triangulation;
};

int main(int argc, const char * argv[]) {

    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::cout << "Reading file : " << input_file << std::endl;

    //// Read file
    std::ifstream input_stream; // stream for input file
    input_stream.open(input_file);
    std::filesystem::path currentPath = std::filesystem::current_path();
    std::cout << "Current working directory: " << currentPath << std::endl;
    if (input_stream.is_open()) {
        std::string line;

        //// Parse line by line
        while (getline(input_stream, line)) {
            std::cout << line << std::endl;

            std::istringstream line_stream(line);
            char line_type; // output the char type in order
            line_stream >> line_type;

            //// Vertex: store all the coordination into vertices
            if (line_type == 'v') {
                vertices.emplace_back(); // append a new vertex
                double x, y, z; // output the double type in order
                line_stream >> x >> y >> z;
                vertices.back().x = x; //.back() returns the last element of the vector = append
                vertices.back().y = y;
                vertices.back().z = z;
            }

            //// Face: store all the vertices into faces
            if (line_type == 'f') {
                faces.emplace_back();
                int v;
                // keep looping until the end of the line
                while (!line_stream.eof()) { // eof = end of file
                    line_stream >> v;
                    faces.back().boundary.emplace_back(v-1);// v-1 because the index starts from 0
                }
            }

        }
    }
    else {
        std::cerr << "Error: unable to open file :" << input_file << std::endl;
        return 1;
    }

////     Print vertices, start with index 0
//  int i = 0;
//  for (auto const &vertex: vertices) {
//    std::cout << "Vertex " << i++ << ": " << "(" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")" << std::endl;
//  }
//
//////    Print faces, start with index 0
//  i = 0;
//  for (auto const &face: faces) {
//    std::cout << "Face " << i++ << ": ";
//    for (auto const &vertex: face.boundary) std::cout << " " << vertex;
//    std::cout << std::endl;
//  }

////     Best fitting planes (to do)
    for (auto &face: faces) { // go through each face
        std::vector<Kernel::Point_3> planes; // create a vector of points
        for (auto const &vertex_id: face.boundary) { // go through each vertex
            planes.emplace_back(vertices[vertex_id].x, vertices[vertex_id].y, vertices[vertex_id].z);
        } // add the vertex to the vector of points
        CGAL::linear_least_squares_fitting_3(planes.begin(),
                                             planes.end(),
                                             face.best_plane,
                                             CGAL::Dimension_tag<0>());// fitting the plane vertices
    }

////     Triangulate faces (to do)
    for (auto &face: faces) {
        std::vector<Kernel::Point_3> projected_boundary;
        for (auto &vertex_id: face.boundary) {
            Kernel::Point_3 p(vertices[vertex_id].x, vertices[vertex_id].y, vertices[vertex_id].z);
            // project the point to the plane
            Kernel::Point_3 projected_p = face.best_plane.projection(p);
            // store the projected point to the vector
            projected_boundary.push_back(projected_p);
            // connect the points in the vector
//            for (int i = 0; i < projected_boundary.size(); ++i) {
//            }

        }
    }
////     Label triangulation (to do)
    for (auto &face: faces) {
        for (auto face_it = face.triangulation.finite_faces_begin(); face_it != face.triangulation.finite_faces_end(); ++face_it) {
            if (face.triangulation.is_infinite(face_it)) {
                face_it->info().interior = false;
            }
            else {
                face_it->info().interior = true;
            }
        }
    }
////     Export triangles (to do)
    std::ofstream output_stream;
    output_stream.open(output_file);
    if (output_stream.is_open()) {
        for (auto const &face: faces) {
            for (auto face_it = face.triangulation.finite_faces_begin(); face_it != face.triangulation.finite_faces_end(); ++face_it) {
                if (face_it->info().interior) {
                    output_stream << "f";
                    for (int i = 0; i < 3; ++i) {
                        Kernel::Point_2 p = face.triangulation.point(face_it, i);
                        output_stream << " " << p.x() << " " << p.y();
                    }
                    output_stream << std::endl;
                }
            }
        }
    }
    else {
        std::cerr << "Error: unable to open file :" << output_file << std::endl;
        return 1;
    }
    return 0;
}
