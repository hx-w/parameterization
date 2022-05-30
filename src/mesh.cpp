#include "mesh.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace glm;

namespace RenderSpace {
    Mesh::~Mesh() {
        _reset();
    }

    Mesh::Mesh(const Mesh& element) {
        _deepcopy(element);
    }

    Mesh& Mesh::operator=(const Mesh& element) {
        _deepcopy(element);
        return *this;
    }

    void Mesh::_reset() {
        vector<Triangle>().swap(m_triangles);
        vector<Vertex>().swap(m_vertices);
    }

    void Mesh::_deepcopy(const Mesh& element) {
        m_triangles.assign(element.m_triangles.begin(), element.m_triangles.end());
        m_vertices.assign(element.m_vertices.begin(), element.m_vertices.end());
    }

    bool Mesh::load_OBJ(const std::string& filename) {
        _reset();
        ifstream ifs(filename);
        if (!ifs.good()) {
            cout << "[ERROR] " << "Can't open file: " << filename << endl;
            return false;
        }

        string line;
        while (getline(ifs, line)) {
            if (line.empty()) {
                continue;
            }
            if (line[0] == '#') {
                continue;
            }
            vector<string> words;
            _split_words(line, words);
            if (words[0] == "v") {
                m_vertices.emplace_back(
                    vec3(stof(words[1]), stof(words[2]), stof(words[3])),
                    vec3(1.0),
                    vec3(0.0)
                );
            }
            else if (words[0] == "f") {
                int v1 = stoi(words[1]) - 1;
                int v2 = stoi(words[2]) - 1;
                int v3 = stoi(words[3]) - 1;
                m_triangles.emplace_back(Triangle(v1, v2, v3));
            }
        }
        ifs.close();

        // 计算法线
        for (int i = 0; i < m_triangles.size(); ++i) {
            Triangle& tri = m_triangles[i];
            vec3 v1 = m_vertices[tri.VertexIdx.x].Position;
            vec3 v2 = m_vertices[tri.VertexIdx.y].Position;
            vec3 v3 = m_vertices[tri.VertexIdx.z].Position;
            vec3 normal = normalize(cross(v2 - v1, v3 - v1));
            m_vertices[tri.VertexIdx.x].Normal += normal;
            m_vertices[tri.VertexIdx.y].Normal += normal;
            m_vertices[tri.VertexIdx.z].Normal += normal;
        }
        cout << "[INFO] obj file loaded: " << filename << " (vt: " << m_vertices.size() << ", tri: " << m_triangles.size() << ")" << endl;
        return true;
    }

    bool Mesh::save_OBJ(const string& filename) {
        ofstream ofs(filename);
        // save obj file by m_vertices and m_triangles
        for (auto& v : m_vertices) {
            ofs << "v " << v.Position.x << " " << v.Position.y << " " << v.Position.z << endl;
        }
        for (auto& tri : m_triangles) {
            ofs << "f " << tri.VertexIdx.x + 1 << " " << tri.VertexIdx.y + 1 << " " << tri.VertexIdx.z + 1 << endl;
        }
        ofs.close();

        cout << "[INFO] obj file saved: " << filename << endl;
        return true;
    }

    void Mesh::_split_words(const string& line, vector<string>& words, const char delim) {
        stringstream ss(line);
        string word;
        while (getline(ss, word, delim)) {
            words.push_back(word);
        }
    }
}
