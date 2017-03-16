#include "Patch.h"

Patch::Patch(string s, bool iris_label)
{
    vector<string> str_list = StringSplit(s, ' ');

    label = Str2Num<int>(str_list[0]);

    if (iris_label) {
        for (size_t i=1;i<str_list.size()-2;i++) {
            info.push_back(Str2Num<int>(str_list[i]));
        }
    }
    else {
        for (size_t i=1;i<4;++i) {
            info.push_back(Str2Num<int>(str_list[i]));
        }
    }

    if (info.size() == 3) shape = Tri;
    else if (info.size() == 4) shape = Rect;
}

ostream& operator<<(ostream& os, Patch &patch)
{
    int node_num = patch.get_shape() == Tri ? 3 : 4;
    for (int i=0;i<node_num;i++) os << patch[i] << " ";
    os << endl;
    os << patch.shape << endl;

    return os;
}
