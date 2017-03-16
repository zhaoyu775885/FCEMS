#include "BasisFunction.h"


BasisFunc::BasisFunc(const Mesh &mesh) :
    g_mesh(mesh), common_edge_num(0)
{
    common_edge_num = g_mesh.get_edge_num();
    for (int i=0;i<common_edge_num;++i) {

        FullBF tmp_bf;
        tmp_bf.type = Rwg;
        tmp_bf.pHalf = new HalfBF;
        tmp_bf.nHalf = new HalfBF;
        const vector<int> tmp_edge = g_mesh.get_edge(i);

        tmp_bf.pHalf->hRwg = new halfrwgBF;
        tmp_bf.pHalf->type = Tri;
        tmp_bf.pHalf->face_label = tmp_edge[2];
        halfrwgBF *tmp_prwg = tmp_bf.pHalf->hRwg;
        tmp_prwg->v0 = g_mesh.get_node(tmp_edge[4]-1);
        tmp_prwg->v1 = g_mesh.get_node(tmp_edge[0]-1);
        tmp_prwg->v2 = g_mesh.get_node(tmp_edge[1]-1);

        Point3D p0(tmp_prwg->v0), p1(tmp_prwg->v1), p2(tmp_prwg->v2);
        tmp_prwg->nvec = get_trigon_nvec(p0, p1, p2);
        tmp_prwg->ctr = get_trigon_ctr(p0, p1, p2);
        tmp_prwg->area = get_trigon_area(p0, p1, p2);
        tmp_prwg->len = (p1-p2).length();
        tmp_prwg->face_label = tmp_edge[2];

        tmp_bf.nHalf->hRwg = new halfrwgBF;
        tmp_bf.nHalf->type = Tri;
        tmp_bf.nHalf->face_label = tmp_edge[3];
        halfrwgBF *tmp_nrwg = tmp_bf.nHalf->hRwg;
        tmp_nrwg->v0 = g_mesh.get_node(tmp_edge[5]-1);
        tmp_nrwg->v1 = g_mesh.get_node(tmp_edge[0]-1);
        tmp_nrwg->v2 = g_mesh.get_node(tmp_edge[1]-1);

        p0 = tmp_nrwg->v0;
        p1 = tmp_nrwg->v1;
        p2 = tmp_nrwg->v2;
        tmp_nrwg->nvec = get_trigon_nvec(p2, p1, p0);
        tmp_nrwg->ctr = get_trigon_ctr(p0, p1, p2);
        tmp_nrwg->area = get_trigon_area(p0, p1, p2);
        tmp_nrwg->len = (p1-p2).length();
        tmp_nrwg->face_label = tmp_edge[3];
//
//        cout << "p: " << tmp_prwg->nvec << endl;
//        cout << tmp_prwg->ctr << endl;
//        cout << "n: " << tmp_nrwg->nvec << endl;
//        cout << tmp_nrwg->ctr << endl;
//        cout << endl;

        bf_database.push_back(tmp_bf);
    }
}
