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
        std::vector<Kernel::Point_2> projected_2d_boundary;
        CGAL::Vector_3<Kernel> normal = face.best_plane.orthogonal_vector(); // normal vector of the plane
        CGAL::Vector_3<Kernel> u = face.best_plane.base1(); // base vector u
        CGAL::Vector_3<Kernel> v = face.best_plane.base2(); // base vector v
        auto first_vertex_id = face.boundary.front(); // first vertex id
        Kernel::Point_3 first_vertex(vertices[first_vertex_id].x, vertices[first_vertex_id].y, vertices[first_vertex_id].z); // first vertex
        Kernel::Point_3 first_projected_vertex = face.best_plane.projection(first_vertex); // first projected vertex
        CGAL::Vector_3<Kernel> first_projected_vector(first_projected_vertex, CGAL::ORIGIN);// first projected vector
        std::cout << "2D plane info: " << u << "," << v << "," << first_projected_vertex << std::endl;
        for (auto &vertex_id: face.boundary) {
            // use the first vertex as the origin
            Kernel::Point_3 p(vertices[vertex_id].x, vertices[vertex_id].y, vertices[vertex_id].z);
            Kernel::Point_3 projected_p = face.best_plane.projection(p);
            CGAL::Vector_3<Kernel> projected_vector(projected_p, CGAL::ORIGIN);
            CGAL::Vector_3<Kernel> difference_vector = projected_vector - first_projected_vector;
            // project the 3D points to 2D
            CGAL::Vector_2<Kernel> projected_2d_vector(difference_vector*u, difference_vector*v);
            Kernel ::Point_2 projected_2d_point(projected_2d_vector.x(), projected_2d_vector.y());
            std::cout << "Projected 2D point: " << projected_2d_point << std::endl;
            projected_2d_boundary.push_back(projected_2d_point);
        }
        // constraint triangulation
        for (int i = 0; i < projected_2d_boundary.size(); ++i) {
            // insert constraint with consecutive vertices
            face.triangulation.insert_constraint(projected_2d_boundary[i],
                                                 projected_2d_boundary[(i+1)%projected_2d_boundary.size()]);
        }
            // check the validity of the triangulation
        assert(face.triangulation.is_valid());
            // triangulate the polygon
        for (auto face_it = face.triangulation.finite_faces_begin();
             face_it != face.triangulation.finite_faces_end(); ++face_it)// go through each face
        {
            if (face.triangulation.is_infinite(face_it)) // check if the face is infinite
            {
                face_it->info().interior = false; // Set the current face's interior flag
            }
            else
            {
                face_it->info().interior = true; // Set the current face's interior flag
            }

        }
        // output the vertex of triangles
        for (auto face_it = face.triangulation.finite_faces_begin();
             face_it != face.triangulation.finite_faces_end(); ++face_it) {
            if (face_it->info().interior) {
                std::cout << "Triangle: ";
                for (int i = 0; i < 3; ++i) {
                    auto vertex = face_it->vertex(i);
                    std::cout << " " << vertex->point();
                }
                std::cout << std::endl;
            }
        }
    }

////     Label triangulation (to do)



////     Export triangles (to do)
    std::ofstream output_stream;
    output_stream.open(output_file);
    // output all original vertices
    for (auto const &vertex: vertices) {
        output_stream << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
    }
}