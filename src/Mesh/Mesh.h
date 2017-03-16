#ifndef MESH_H
#define MESH_H

#include "../settings.h"
#include "../Basics/routine.h"
#include "../Math/Point3d.h"
#include "Patch.h"

class Mesh
{
public:
    Mesh(const string &iris_mesh);
    Mesh(const string &feko_node, const string &feko_edge, const string &feko_face);

    const Point3D &get_node(int n) const { return node_info[n]; }
    const Patch   &get_face(int n) const { return face_info[n]; }
    const vector<int> &get_edge(int n) const { return edge_info[n]; }

    int get_node_num() const {return node_num;}
    int get_edge_num() const {return edge_num;}
    int get_face_num() const {return face_num;}

    void readNode_feko(const string &filename);
    void readEdge_feko(const string &filename);
    void readFace_feko(const string &filename);

    void display_mesh_info() const;

private:
    int node_num;
    int edge_num;
    int face_num;

    vector<Point3D> node_info;
    vector<Patch>   face_info;
    vector< vector<int> > edge_info;
};





#endif // MESH_H
