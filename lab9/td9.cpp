#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "classes.h"
#include <set>

#define M_PI 3.14159265358979323846

class TriangleIndex
{
   public:
    TriangleIndex(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1,
                  int nj = -1, int nk = -1, int uvi = -1, int uvj = -1,
                  int uvk = -1, int group = -1, bool added = false)
        : vtxi(vtxi),
          vtxj(vtxj),
          vtxk(vtxk),
          uvi(uvi),
          uvj(uvj),
          uvk(uvk),
          ni(ni),
          nj(nj),
          nk(nk),
          group(group){};
    int vtxi, vtxj, vtxk;  // indices within the vertex coordinates array
    int uvi, uvj, uvk;     // indices within the uv coordinates array
    int ni, nj, nk;        // indices within the normals array
    int group;             // face group
};

class Edge
{
   public:
    Edge(int a = 0, int b = 0)
    {
        if (a <= b)
        {
            vtxA = a;
            vtxB = b;
        }
        else
        {
            vtxA = b;
            vtxB = a;
        }
    }

    bool operator==(const Edge &other)
    {
        return (vtxA == other.vtxA && vtxB == other.vtxB);
    }

    int vtxA;
    int vtxB;
};

bool operator<(const Edge &current, const Edge &other)
{
    if (current.vtxA < other.vtxA)
    {
        return true;
    }
    if (current.vtxA > other.vtxA)
    {
        return false;
    }
    return current.vtxB < other.vtxB;
}

class TriangleMesh
{
   public:
    TriangleMesh() {}
    ~TriangleMesh() {}

    void readOBJ(const char *obj)
    {
        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0],
                           &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndex t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0,
                            &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0,
                                &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2,
                                    &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n",
                                        &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3,
                                &offset);
                    TriangleIndex t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3,
                                        &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
    }

    void writeOBJ(const char *obj)
    {
        FILE *f = fopen(obj, "w+");

        for (int i = 0; i < vertices.size(); i++)
        {
            fprintf(f, "v %f %f %f\n", vertices[i][0], vertices[i][1],
                    vertices[i][2]);
        }

        for (int i = 0; i < indices.size(); i++)
        {
            fprintf(f, "f %u %u %u\n", indices[i].vtxi + 1, indices[i].vtxj + 1,
                    indices[i].vtxk + 1);
        }

        fclose(f);
    }

    //PG 109!!!!

      std::map<Edge, std::vector<int>> buildEdgesToTriangles() {
    std::map<Edge, std::vector<int>> edges_to_triangles;
    int n_triangles = indices.size();
    for (int i = 0; i < n_triangles; ++i) {
        edges_to_triangles[Edge(indices[i].vtxi, indices[i].vtxj)].push_back(i);
        edges_to_triangles[Edge(indices[i].vtxj, indices[i].vtxk)].push_back(i);
        edges_to_triangles[Edge(indices[i].vtxk, indices[i].vtxi)].push_back(i);
    }
    return edges_to_triangles;
    }

    std::map<int, std::set<int>> buildAdjacencyList() {
        std::map<int, std::set<int>> adjacency_list;
        int n_triangles = indices.size();
        for (int i = 0; i < n_triangles; ++i) {
            adjacency_list[indices[i].vtxi].insert(indices[i].vtxj);
            adjacency_list[indices[i].vtxi].insert(indices[i].vtxk);
            adjacency_list[indices[i].vtxj].insert(indices[i].vtxi);
            adjacency_list[indices[i].vtxj].insert(indices[i].vtxk);
            adjacency_list[indices[i].vtxk].insert(indices[i].vtxi);
            adjacency_list[indices[i].vtxk].insert(indices[i].vtxj);
        }
        return adjacency_list;
    }

    std::vector<Edge> findBoundaryEdges(const std::map<Edge, std::vector<int>>& edges_to_triangles) {
        std::vector<Edge> boundary_edges;
        for (const auto& [edge, triangles] : edges_to_triangles) {
            if (triangles.size() == 1) {
                boundary_edges.push_back(edge);
            }
        }
        return boundary_edges;
    }

    std::vector<int> orderBoundaryVertices(const std::vector<Edge>& boundary_edges) {
        std::vector<int> ordered_boundary_vertices;
        int v_first_vertex = boundary_edges[0].vtxA;
        int second_last_vertex = v_first_vertex;
        int last_vertex = boundary_edges[0].vtxB;

        ordered_boundary_vertices.push_back(v_first_vertex);
        ordered_boundary_vertices.push_back(last_vertex);

        while (last_vertex != v_first_vertex) {
            auto it = std::find_if(boundary_edges.begin(), boundary_edges.end(), [&](const Edge& edge) {
                return (edge.vtxA == last_vertex && edge.vtxB != second_last_vertex) ||
                    (edge.vtxB == last_vertex && edge.vtxA != second_last_vertex);
            });
            if (it == boundary_edges.end()) {
                throw std::runtime_error("No edge found");
            }
            if (it->vtxA == last_vertex) {
                second_last_vertex = last_vertex;
                last_vertex = it->vtxB;
                ordered_boundary_vertices.push_back(it->vtxB);
            } else {
                second_last_vertex = last_vertex;
                last_vertex = it->vtxA;
                ordered_boundary_vertices.push_back(it->vtxA);
            }
        }

        return ordered_boundary_vertices;
    }

    double computeBoundaryLength(const std::vector<int>& ordered_boundary_vertices) {
        double s = 0.0;
        for (size_t i = 0; i < ordered_boundary_vertices.size() - 1; ++i) {
            s += (vertices[ordered_boundary_vertices[i + 1]] - vertices[ordered_boundary_vertices[i]]).norm();
        }
        return s;
    }

    std::vector<Vector> layoutBoundaryVerticesOnCircle(const std::vector<int>& ordered_boundary_vertices, double s) {
        double cs = 0.0;
        std::vector<Vector> v_prime(vertices.size());
        for (size_t i = 0; i < ordered_boundary_vertices.size() - 1; ++i) {
            double theta_i = 2.0 * M_PI * cs / s;
            v_prime[ordered_boundary_vertices[i]] = Vector(cos(theta_i), sin(theta_i), 0.0);
            cs += (vertices[ordered_boundary_vertices[i + 1]] - vertices[ordered_boundary_vertices[i]]).norm();
        }
        return v_prime;
    }

    std::vector<Vector> layoutInternalVertices(std::vector<Vector>& v_prime, const std::map<int, std::set<int>>& adjacency_list, const std::vector<int>& ordered_boundary_vertices, int n_iter) {
        for (int iter = 0; iter < n_iter; ++iter) {
            std::cout << "Iteration " << iter << std::endl;
            std::vector<Vector> v_next(v_prime);
            for (size_t i = 0; i < vertices.size(); ++i) {
                if (std::find(ordered_boundary_vertices.begin(), ordered_boundary_vertices.end(), i) != ordered_boundary_vertices.end()) {
                    continue;
                }
                Vector sum(0.0, 0.0, 0.0);
                int count = 0;
                for (const auto& neighbor : adjacency_list.at(i)) { // Use 'at' instead of '[]' to access the set of neighbors
                    sum = sum + v_prime[neighbor];
                    count++;
                }
                v_next[i] = sum / count;
            }
            v_prime = v_next;
        }
        return v_prime;
    }

    // Milos Oundjian helped me grately with this
    TriangleMesh tutte(int n_iter = 1000) {
        std::map<Edge, std::vector<int>> edges_to_triangles = buildEdgesToTriangles();
        std::map<int, std::set<int>> adjacency_list = buildAdjacencyList();
        std::vector<Edge> boundary_edges = findBoundaryEdges(edges_to_triangles);
        std::vector<int> ordered_boundary_vertices = orderBoundaryVertices(boundary_edges);
        double s = computeBoundaryLength(ordered_boundary_vertices);
        std::vector<Vector> v_prime = layoutBoundaryVerticesOnCircle(ordered_boundary_vertices, s);
        v_prime = layoutInternalVertices(v_prime, adjacency_list, ordered_boundary_vertices, n_iter);
        vertices = v_prime;
        return *this;
    }

    

    std::vector<TriangleIndex> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
};

int main()
{
    TriangleMesh m;
    m.readOBJ("mask.obj");
    TriangleMesh result = m.tutte(1000);
    result.writeOBJ("ouput.obj");

    return 0;
}