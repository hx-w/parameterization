﻿#include "parameterization.h"
#include <omp.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <utility>
#include <chrono>
#include <ctime>

using namespace std;
using namespace glm;

namespace RenderSpace {
Parameterization::Parameterization(Mesh* uns_mesh, Mesh* param_mesh, Mesh* str_mesh)
    : m_uns_mesh(uns_mesh), m_param_mesh(param_mesh), m_str_mesh(str_mesh), m_scale(10) {}

Parameterization::~Parameterization() {}

void Parameterization::parameterize() {
    vector<OrderedEdge> edge_bound;
    vector<OrderedEdge> edge_inner;
    // 先获取边缘和非边缘边
    _remark_edges(edge_bound, edge_inner);
    // 需要将边缘边有序遍历 v1->v3, v3->v2, v2->v1
    _topology_reorder(edge_bound);
    // 参数平面 边缘点 根据edge_bound 顺序计算得到
    vector<vec2> param_bound;
    _parameterize_bound(edge_bound, param_bound);
    cout << "param_bound mapping finished" << endl;
    // 初始化weights
    // weight[(i, j)] = weight[(j, i)], 故存储时令i < j
    // 这里元素量太大，可以使用map进行优化，但是map对于key的查找有问题(?)
    _init_weights(edge_bound, edge_inner);
    // 将边集 转换为 点集 保存映射的顺序关系
    // 其中边缘点的顺序不能变，应与param_bound一致
    vector<int> vt_bound;
    vector<int> vt_inner;
    _convert_edge_to_vertex(move(edge_bound), move(edge_inner), vt_bound,
                            vt_inner);
    // 解方程组
    vector<vec2> param_inner;
    cout << "solving Laplacian equation" << endl;
    _solve_Laplacian_equation(vt_inner, vt_inner, param_inner, vt_inner,
                              vt_bound, param_bound);
    cout << "building mesh" << endl;
    _build_param_mesh(vt_inner, vt_bound, param_inner, param_bound);
}

void Parameterization::resample(uint32_t num_samples) {
    Mesh sample_mesh;
    auto& uns_vertices = m_uns_mesh->get_vertices();
    auto& param_vertices = m_param_mesh->get_vertices();
    auto& sample_vertices = sample_mesh.get_vertices();
    auto& str_vertices = m_str_mesh->get_vertices();
    auto& str_trias = m_str_mesh->get_triangles();

    for (auto ir = 0; ir < num_samples; ++ir) {
        for (auto ic = 0; ic < num_samples; ++ic) {
            // 在参数平面上的点
            auto x = m_scale * (ic + 0.5) / num_samples - m_scale / 2;
            auto y = m_scale * (ir + 0.5) / num_samples - m_scale / 2;
            sample_vertices.push_back(Vertex(vec3(x, y, 0), vec3(1.0), vec3(0.0)));
            // 逆映射到三维空间
            const Triangle spot_trias = _which_trias_in(vec2(x, y));
            auto tot_area = _trias_area(param_vertices[spot_trias.VertexIdx.x].Position,
                                        param_vertices[spot_trias.VertexIdx.y].Position,
                                        param_vertices[spot_trias.VertexIdx.z].Position);
            auto vi_area = _trias_area(vec3(x, y, 0),
                                       param_vertices[spot_trias.VertexIdx.y].Position,
                                       param_vertices[spot_trias.VertexIdx.z].Position);
            auto vj_area = _trias_area(param_vertices[spot_trias.VertexIdx.x].Position,
                                       vec3(x, y, 0),
                                       param_vertices[spot_trias.VertexIdx.z].Position);
            auto vk_area = _trias_area(param_vertices[spot_trias.VertexIdx.x].Position,
                                        param_vertices[spot_trias.VertexIdx.y].Position,
                                        vec3(x, y, 0));
            vec3 str_point = (vi_area * uns_vertices[spot_trias.VertexIdx.x].Position +
                              vj_area * uns_vertices[spot_trias.VertexIdx.y].Position +
                              vk_area * uns_vertices[spot_trias.VertexIdx.z].Position) / tot_area;
            
            // shift
            str_point += vec3(10.0, 0.0, 0.0);
            str_vertices.push_back(Vertex(str_point, vec3(1.0), vec3(0.0)));
            /**
             *  retopology
             *  [idx-max_col-1] ----- [idx-max_col]
             *   |            /        |
             *  [idx-1] ------------- [idx] 
             */
            if (ir > 0 && ic > 0) {
                auto idx = sample_vertices.size() - 1;
                str_trias.push_back(Triangle(idx, idx - num_samples, idx - 1));
                str_trias.push_back(Triangle(idx - 1, idx - num_samples, idx - num_samples - 1));
            }
        }
    }

    sample_mesh.save_OBJ("models/sample.obj");
}

void Parameterization::_remark_edges(vector<OrderedEdge>& edge_bound,
                                     vector<OrderedEdge>& edge_inner) {
    set<OrderedEdge> edge_bound_set;
    set<OrderedEdge> edge_inner_set;
    // 构造Edge集合，判断只与一个三角形相邻的边
    unordered_map<OrderedEdge, int, pair_hash> edge_count_map;
    const auto& trias = m_uns_mesh->get_triangles();
    for (auto& tri : trias) {
        for (int i = 0; i < 3; ++i) {
            int vidx = std::min(tri.VertexIdx[i], tri.VertexIdx[(i + 1) % 3]);
            int vidx_next = std::max(tri.VertexIdx[i], tri.VertexIdx[(i + 1) % 3]);
            OrderedEdge edge(vidx, vidx_next);
            if (edge_count_map.find(edge) == edge_count_map.end()) {
                edge_count_map[edge] = 1;
            } else {
                edge_count_map[edge]++;
            }
        }
    }
    for (auto& [edge, count] : edge_count_map) {
        if (count == 1) {
            edge_bound_set.insert(edge);
        } else if (count == 2) {
            edge_inner_set.insert(edge);
        } else {
            cout << "Error: edge count is not 1 or 2" << endl;
        }
    }
    edge_bound.assign(edge_bound_set.begin(), edge_bound_set.end());
    edge_inner.assign(edge_inner_set.begin(), edge_inner_set.end());
}

void Parameterization::_topology_reorder(vector<OrderedEdge>& edge_bound) {
    if (edge_bound.empty())
        return;
    // 对边缘边进行 拓扑关系排序
    // 先构造一个邻接表
    unordered_map<int, vector<int>> adj_list;
    for (auto& edge : edge_bound) {
        adj_list[edge.first].push_back(edge.second);
        adj_list[edge.second].push_back(edge.first);
    }

    // 理论上， 所有边缘点构成一个环，所有adj个数都为2
    vector<OrderedEdge> edge_bound_reorder;
    // 从任一点触发，构造一个拓扑序列
    int vidx = edge_bound[0].first;
    int vidx_reserved = vidx;
    set<int> visited;

    // 顺带计算边缘总长度
    const auto& vertices = m_uns_mesh->get_vertices();
    m_bound_length = 0;

    while (visited.count(vidx) == 0) {
        visited.insert(vidx);
        int vidx_next = adj_list[vidx][0];
        if (visited.count(adj_list[vidx][0]) != 0) {
            // 判断是否首尾相接
            if (visited.count(adj_list[vidx][1]) != 0) {
                vidx_next = vidx_reserved;
            } else {
                vidx_next = adj_list[vidx][1];
            }
        }
        edge_bound_reorder.push_back(OrderedEdge(vidx, vidx_next));
        m_bound_length +=
            length(vertices[vidx].Position - vertices[vidx_next].Position);
        vidx = vidx_next;
    }
    edge_bound.swap(edge_bound_reorder);
}

void Parameterization::_parameterize_bound(vector<OrderedEdge>& edge_bound,
                                           vector<vec2>& param_bound) {
    if (edge_bound.empty())
        return;
    // 参数平面 边缘点 根据edge_bound 顺序计算得到
    const auto& vertices = m_uns_mesh->get_vertices();
    // 三维空间中网格的边缘会被映射到二维参数平面的单位圆/或正方形边缘
    // param_x^j = sin(\theta^j)
    // param_y^j = cos(\theta^j)
    // \theta^j = 2 * \pi (\sum_{i=1}^{j} (vb_{i + 1} - vb_{i})) /
    // m_bound_length 其中 vb_{i + 1} - vb_{i} 为边缘点的距离
    param_bound.clear();
    float _accumulate_length = 0.0;
    bool disturbed_1 = false;
    bool disturbed_2 = false;
    bool disturbed_3 = false;
    for (auto& edge : edge_bound) {
        _accumulate_length +=
            length(vertices[edge.first].Position - vertices[edge.second].Position);
        // float _theta = 2 * M_PI * _accumulate_length / m_bound_length;
        // param_bound.push_back(vec2(sin(_theta) * 10, cos(_theta) * 10));

        /**
         * mapping to rectangle bound
         * make sure every corner of rect is mapped to a vertex
         */
        vec2 bound_point(0.0, 0.0);
        float ratio = _accumulate_length / m_bound_length;
        if (ratio < 0.25) {
            bound_point.x = -(m_scale / 2) + m_scale * (ratio / 0.25);
            bound_point.y = -(m_scale / 2);
        } else if (ratio < 0.5) {
            if (!disturbed_1) {
                disturbed_1 = true;
                ratio = 0.25;
            }
            bound_point.x = (m_scale / 2);
            bound_point.y = -(m_scale / 2) + m_scale * ((ratio - 0.25) / 0.25);
        } else if (ratio < 0.75) {
            if (!disturbed_2) {
                disturbed_2 = true;
                ratio = 0.5;
            }
            bound_point.x = (m_scale / 2) - m_scale * ((ratio - 0.5) / 0.25);
            bound_point.y = (m_scale / 2);
        } else {
            if (!disturbed_3) {
                disturbed_3 = true;
                ratio = 0.75;
            }
            bound_point.x = -(m_scale / 2);
            bound_point.y = (m_scale / 2) - m_scale * ((ratio - 0.75) / 0.25);
        }
        param_bound.push_back(bound_point);
    }
}

void Parameterization::_init_weights(const vector<OrderedEdge>& edge_bound,
                                     const vector<OrderedEdge>& edge_inner) {
    /**
     * weights 的计算满足以下约束
     * 1. weights的key记作(vi, vj)，满足vi <= vj
     * 2. vi与vj不会都是边缘点
     *
     * weights[(vi, vj)] 的计算方法如下
     * 存在以边(vi, vj)为公共边的两个三角形△(vi, vl, vj)，△(vj, vk, vi)
     * 定义角度\alpha(vi, vj) = angle(vec(vi, vl), vec(vl, vj))
     *        \alpha(vj, vi) = angle(vec(vj, vk), vec(vk, vi))
     * weights[(vi, vj)] = [cot(\alpha(vi, vj)) + cot(\alpha(vj, vi))] / 2
     */
    const auto& vertices = m_uns_mesh->get_vertices();
    const auto& triangles = m_uns_mesh->get_triangles();
    unordered_set<Triangle, trias_hash> trias_set(triangles.begin(),
                                                  triangles.end());
    /** 存在这样一个事实：
     *  vi的邻接点中有vk, vk的邻接点中有vj，可能并不存在三角形(vi, vj, vk)
     * (vi) --------------(vk)
     *  | \           /  /
     *  |  \        /  /
     *  |    -(vl)-  /
     *  |     /    /
     *  |   /    /
     *  | /   /
     * (vj)-
     */
    // 构造邻接表 快速索引
    unordered_map<int, set<int>> adj_list;
    vector<OrderedEdge> tot_edge;
    tot_edge.insert(tot_edge.begin(), edge_bound.begin(), edge_bound.end());
    tot_edge.insert(tot_edge.begin(), edge_inner.begin(), edge_inner.end());
    for (auto& edge : tot_edge) {
        adj_list[edge.first].insert(edge.second);
        adj_list[edge.second].insert(edge.first);
    }

    // weights只需从inner构造 [X] bound 也需要
    for (auto edge : tot_edge) {
        int vi = std::min(edge.first, edge.second);
        int vj = std::max(edge.first, edge.second);
        vector<int> adj_vt;
        for (auto& adj_v : adj_list[vi]) {
            // 判断(adj_v, vj, vi)是否构成三角形
            if (adj_v != vj && trias_set.count(Triangle(vi, vj, adj_v)) != 0) {
                adj_vt.push_back(adj_v);
            }
        }
        if (adj_vt.size() > 2) {
            cout << "Error: edge_inner " << vi << " " << vj << ": "
                 << adj_vt.size() << endl;
            continue;
        }

        float _weight = 0.0f;
        for (auto vk : adj_vt) {
            _weight += _cot(_angle_between(vertices[vi].Position,
                                           vertices[vj].Position,
                                           vertices[vk].Position));
        }

        _weight /= adj_vt.size();

        if (fabs(_weight) < 1e-6) {
            cout << "too low: edge_inner " << vi << " " << vj << ": " << _weight
                 << endl;
        }
        m_weights[OrderedEdge(vi, vj)] = _weight;
        // weights中 i=j无意义，但是可以预存ij相等的情况，方便Laplacian
        // matrix的计算 默认值是0
        m_weights_diag[vi] += _weight;
        m_weights_diag[vj] += _weight;
    }
}

void Parameterization::_convert_edge_to_vertex(vector<OrderedEdge>&& edge_bound,
                                               vector<OrderedEdge>&& edge_inner,
                                               vector<int>& vt_bound,
                                               vector<int>& vt_inner) {
    if (edge_bound.empty() || edge_inner.empty()) {
        return;
    }
    vt_bound.clear();
    vt_inner.clear();
    // 计算vt_bound
    // 由于边缘边拓扑有序，所以尽量避免在vt_bound中查重操作
    const auto sz = edge_bound.size();
    // 最后一个边包含的两个点都应已经被包含在vt_bound中
    vt_bound.emplace_back(edge_bound[0].first);
    for (size_t eidx = 0; eidx < sz - 1; ++eidx) {
        const int _last_vt = vt_bound.back();
        if (edge_bound[eidx].first == _last_vt) {
            vt_bound.emplace_back(edge_bound[eidx].second);
        } else if (edge_bound[eidx].second == _last_vt) {
            vt_bound.emplace_back(edge_bound[eidx].first);
        } else {
            cout << "Error: edge_bound is not topology sorted" << endl;
        }
    }
    // 计算vt_inner
    // 先将vt_bound 构造成集合，查找效率高
    // (vt_bound元素少，使用基于hash的unordered_set)
    unordered_set<int> vt_bound_set(vt_bound.begin(), vt_bound.end());
    set<int> vt_inner_set;
    for (auto& edge : edge_inner) {
        if (vt_bound_set.count(edge.first) == 0) {
            vt_inner_set.insert(edge.first);
        }
        if (vt_bound_set.count(edge.second) == 0) {
            vt_inner_set.insert(edge.second);
        }
    }

    vt_inner.assign(vt_inner_set.begin(), vt_inner_set.end());
    // 释放内存
    vector<OrderedEdge>().swap(edge_bound);
    vector<OrderedEdge>().swap(edge_inner);
}

void Parameterization::_solve_Laplacian_equation(
    const vector<int>& r_idx_1,
    const vector<int>& c_idx_1,
    vector<vec2>& f_1,  // 结果保存在这里
    const vector<int>& r_idx_2,
    const vector<int>& c_idx_2,
    const vector<vec2>& f_2) {
    // 约束 c_idx_2.size() == f_2.size()
    //     r_idx_1.size() == c_idx_1.size() 即矩阵1为方阵
    // 变量初始化
    const int mat1_row_count = r_idx_1.size();
    const int mat1_col_count = c_idx_1.size();
    const int mat2_row_count = r_idx_2.size();
    const int mat2_col_count = c_idx_2.size();
    assert(mat2_col_count == f_2.size());
    // L(r1, c1) * f1 = -L(r2, c2) * f2
    // 令 _value_mat = -L(r2, c2) * f2

    vector<vec2> _value_mat;
    for (int ir = 0; ir < mat2_row_count; ++ir) {
        vec2 _row_vec(0.0, 0.0);
        for (int ic = 0; ic < mat2_col_count; ++ic) {
            float _v = _Laplacian_val(r_idx_2[ir], c_idx_2[ic]);
            _row_vec += _v * f_2[ic];
        }
        _value_mat.push_back(-_row_vec);
    }
    // 设置迭代初值 (0.0, 0.0)
    f_1.resize(mat1_col_count, vec2(0.5f, 0.5f));
    // 进行迭代求解

    Jacobi_Iteration(r_idx_1, c_idx_1, f_1, _value_mat, 0.0001f);
}

void Parameterization::Jacobi_Iteration(const vector<int>& r_idx,
                                              const vector<int>& c_idx,
                                              vector<vec2>& f,
                                              const vector<vec2>& b,
                                              const float epsilon  // 允许的误差
) {
    const int row_max = r_idx.size();
    const int col_max = c_idx.size();
    const int f_max = f.size();
    // row_max == col_max == f_max
    assert(row_max == col_max);
    assert(row_max == f_max);

    const int _max_iter = 100;  // 最大迭代次数
    for (int _iter_count = 0; _iter_count < _max_iter; ++_iter_count) {
        float _residual = 0.0f;
        auto start = chrono::system_clock::now();
        vector<vec2> _new_f(f);
        // 对于x_{f_max}^{_iter_count}
#pragma omp parallel for reduction(+:_residual)
        for (int ir = 0; ir < f_max; ++ir) {
            vec2 _val(0.0, 0.0);
            for (int ic = 0; ic < f_max; ++ic) {
                if (ir == ic)
                    continue;
                float _lp = _Laplacian_val(r_idx[ir], c_idx[ic]);
                _val += _lp * f[ic];
            }
            float _iv = -1.0 / _Laplacian_val(r_idx[ir], c_idx[ir]);
            _val = (_val - b[ir]) * _iv;
            _residual += sqrt(length(f[ir] - _val));
            _new_f[ir] = _val;
        }
        f.assign(_new_f.begin(), _new_f.end());
        auto end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        time_t end_time = chrono::system_clock::to_time_t(end);
        cout << ">> " << ctime(&end_time) << _iter_count << "/" << _max_iter << " ==> " << _residual << "  | cost " << elapsed_seconds.count() << endl;
        if (_residual < epsilon) {
            break;
        }
    }
}

void Parameterization::_build_param_mesh(const vector<int>& vt_inner,
                                   const vector<int>& vt_bound,
                                   const vector<vec2>& param_inner,
                                   const vector<vec2>& param_bound) {
    cout << "[DEBUG] rebuilding" << endl;
    // 对vt_inner, vt_bound构建倒排索引
    map<int, int> vt_inner_idx;
    map<int, int> vt_bound_idx;
    for (int i = 0; i < vt_inner.size(); ++i) {
        vt_inner_idx[vt_inner[i]] = i;
    }
    for (int i = 0; i < vt_bound.size(); ++i) {
        vt_bound_idx[vt_bound[i]] = i;
    }

    auto& tar_vertices = m_param_mesh->get_vertices();
    auto& tar_tris = m_param_mesh->get_triangles();
    const auto& ori_vertices = m_uns_mesh->get_vertices();
    const auto& ori_tris = m_uns_mesh->get_triangles();
    tar_vertices.clear();

    // 三角面片索引应相同
    tar_tris.assign(ori_tris.begin(), ori_tris.end());

    // 重设顶点位置
    const int sz = ori_vertices.size();
    for (int i = 0; i < sz; ++i) {
        vec2 _v(0.0, 0.0);
        if (vt_inner_idx.find(i) != vt_inner_idx.end()) {
            _v = param_inner[vt_inner_idx[i]];
        } else if (vt_bound_idx.find(i) != vt_bound_idx.end()) {
            _v = param_bound[vt_bound_idx[i]];
        } else {
            cout << "[ERROR] 发现非法顶点索引" << endl;
        }
        tar_vertices.push_back(
            Vertex(vec3(_v.x, _v.y, 0.0), vec3(1.0), vec3(1.0)));
    }
}

float Parameterization::_Laplacian_val(int i, int j) {
    if (i > j)
        swap(i, j);
    if (i != j) {
        auto iter = m_weights.find(OrderedEdge(i, j));
        if (iter == m_weights.end()) {
            return 0.0f;
        } else {
            return -iter->second;
        }
    } else {
        return m_weights_diag.find(i)->second;
    }
    return 0;
}

float Parameterization::_cot(float rad) const {
    return cos(rad) / sin(rad);
}

float Parameterization::_angle_between(const vec3& va,
                                       const vec3& vb,
                                       const vec3& ori) const {
    // calculate angle between va-ori and vb-ori by glm
    vec3 _va = normalize(va - ori);
    vec3 _vb = normalize(vb - ori);
    float _cos = dot(_va, _vb);
    float _rad = acos(_cos);
    return _rad;
}

float Parameterization::_trias_area(const vec3& v0, const vec3& v1, const vec3& v2) const {
    float _a = length(v1 - v0);
    float _b = length(v2 - v0);
    float _c = length(v2 - v1);
    float _s = (_a + _b + _c) / 2.0f;
    return sqrt(_s * (_s - _a) * (_s - _b) * (_s - _c));
}

const Triangle Parameterization::_which_trias_in(const vec2& pos) const {
    const auto& _vertices = m_param_mesh->get_vertices();
    const auto& _tris = m_param_mesh->get_triangles();
    const int _sz = _vertices.size();
    for (const auto& tri : _tris) {
        const vec3& _v0 = _vertices[tri.VertexIdx.x].Position;
        const vec3& _v1 = _vertices[tri.VertexIdx.y].Position;
        const vec3& _v2 = _vertices[tri.VertexIdx.z].Position;
        // judge whether pos in triangle [_v0, _v1, _v2]
        auto s = (_v0.x - _v2.x) * (pos.y - _v2.y) - (_v0.y - _v2.y) * (pos.x - _v2.x);
        auto t = (_v1.x - _v0.x) * (pos.y - _v0.y) - (_v1.y - _v0.y) * (pos.x - _v0.x);

        if ((s < 0) != (t < 0) && s != 0 && t != 0) {
            continue;
        }

        auto d = (_v2.x - _v1.x) * (pos.y - _v1.y) - (_v2.y - _v1.y) * (pos.x - _v1.x);
        if (d == 0 || (d < 0) == (s + t <= 0)) {
            return tri;
        }
    }
    return Triangle(0.0, 0.0, 0.0);
}

}  // namespace RenderSpace
