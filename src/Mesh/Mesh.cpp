#include "Mesh.h"

Mesh::Mesh(const string &feko_node, const string &feko_edge, const string &feko_face)
{
    readNode_feko(feko_node);
    readEdge_feko(feko_edge);
    readFace_feko(feko_face);
}

void
Mesh::readNode_feko(const string &filename)
{
    string line_str;
    ifstream fp(filename.c_str());
    if (!fp) {
        cerr << "error open " << filename << "... (in node info input process )" << endl;
        exit(1);
    }
    else {
        getline(fp, line_str);
        vector<string> vec_str = StringSplit(line_str, ' ');
        node_num = Str2Num<int>( vec_str[0] );
        for (int j=0;j<node_num;j++) {
            getline(fp, line_str);
            Point3D point(line_str, false);
//            cout << point << endl;
            node_info.push_back(point);
        }
    }
    fp.close();
}

void
Mesh::readEdge_feko(const string &filename)
{
    string line_str;
    ifstream fp(filename.c_str());
    if (!fp) {
        cerr << "error open " << filename << "... (in face info input process )" << endl;
        exit(1);
    }
    else {
        getline(fp, line_str);
        vector<string> vec_str = StringSplit(line_str, ' ');

        edge_num = Str2Num<int>( vec_str[0] );
        int item_num = 6;
        for (int i=0;i<edge_num;++i) {
            getline(fp, line_str);
            vec_str = StringSplit(line_str, ' ');
            vector<int> tmp_edge;
            for (int j=0;j<item_num;++j) tmp_edge.push_back(Str2Num<int>( vec_str[j] ));
            edge_info.push_back(tmp_edge);
        }
    }

    fp.close();
}

void
Mesh::readFace_feko(const string &filename)
{
    string line_str;
    ifstream fp(filename.c_str());
    if (!fp) {
        cerr << "error open " << filename << "... (in face info input process )" << endl;
        exit(1);
    }
    else {
        getline(fp, line_str);
        vector<string> vec_str = StringSplit(line_str, ' ');
        face_num = Str2Num<int>( vec_str[0] );
        for (int i=0;i<face_num;++i) {
            getline(fp, line_str);
            Patch patch(line_str, false);
            face_info.push_back(patch);
        }
    }

    fp.close();
}

void
Mesh::display_mesh_info() const
{
    cout << "node number: " << get_node_num() << endl;
    cout << "edge number: " << get_edge_num() << endl;
    cout << "face number: " << get_face_num() << endl;
    cout << "***************************************" << endl;
}
