#ifndef Mesh_H
#define Mesh_H
#include <mutex>
#include <vector>
#include <string>
#include <algorithm>
#include "glm/glm.hpp"

namespace RenderSpace {
    typedef std::pair<int, int> OrderedEdge; // v1, v2; v1<v2

    struct Vertex {
        Vertex() = default;
        Vertex(const glm::vec3& pos, const glm::vec3& clr, const glm::vec3& nml) :
            Position(pos), Color(clr), Normal(nml) { }
        glm::vec3 Position;
        glm::vec3 Color;
        glm::vec3 Normal;
    };

    struct Triangle {
        Triangle() = default;
        // 顶点顺序不能变
        Triangle(int v0, int v1, int v2):
            VertexIdx(glm::uvec3(v0, v1, v2)) { }
        bool operator==(const Triangle& other) const {
            return (VertexIdx.x == other.VertexIdx.x && VertexIdx.y == other.VertexIdx.y && VertexIdx.z == other.VertexIdx.z) \
                || (VertexIdx.x == other.VertexIdx.x && VertexIdx.y == other.VertexIdx.z && VertexIdx.z == other.VertexIdx.y) \
                || (VertexIdx.x == other.VertexIdx.y && VertexIdx.y == other.VertexIdx.z && VertexIdx.z == other.VertexIdx.x) \
                || (VertexIdx.x == other.VertexIdx.y && VertexIdx.y == other.VertexIdx.x && VertexIdx.z == other.VertexIdx.z) \
                || (VertexIdx.x == other.VertexIdx.z && VertexIdx.y == other.VertexIdx.x && VertexIdx.z == other.VertexIdx.y) \
                || (VertexIdx.x == other.VertexIdx.z && VertexIdx.y == other.VertexIdx.y && VertexIdx.z == other.VertexIdx.x);
        }
        bool operator<(const Triangle& other) const {
            return VertexIdx.x < other.VertexIdx.x;
        }

        int vt_in(int vt) const {
            if (VertexIdx.x == vt) return 0;
            if (VertexIdx.y == vt) return 1;
            if (VertexIdx.z == vt) return 2;
            return -1;
        }

        bool is_neighbor(const Triangle& tri) {
            int x_ret = vt_in(tri.VertexIdx.x);
            int y_ret = vt_in(tri.VertexIdx.y);
            int z_ret = vt_in(tri.VertexIdx.z);
            return ((x_ret >= 0 && y_ret >= 0) || (x_ret >= 0 && z_ret >= 0) || (y_ret >= 0 && z_ret >= 0));
        }

        glm::uvec3 VertexIdx;
    };

    // base class for Mesh objects
    class Mesh {
    public:
        Mesh() = default;
        ~Mesh();

        Mesh(const Mesh&);
        Mesh& operator=(const Mesh&);

        bool load_OBJ(const std::string& filename);
        bool save_OBJ(const std::string& filename);

        std::vector<Triangle>& get_triangles() {
            return m_triangles;
        }
        std::vector<Vertex>& get_vertices() {
            return m_vertices;
        }

    private:
        void _deepcopy(const Mesh& element);
        void _split_words(const std::string& line, std::vector<std::string>& words, const char delim=' ');
        void _reset();

    private:
        std::vector<Triangle> m_triangles;
        std::vector<Vertex> m_vertices;
    };
}

#endif
