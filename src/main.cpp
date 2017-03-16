#include <iostream>
#include "Mesh/Mesh.h"
#include "BasisFunction.h"
#include "BasicBEM.h"
#include "PostProcess.h"

#include "h2lib/cluster.h"
#include "h2lib/block.h"
#include "h2lib/hmatrix.h"
#include "h2lib/harith.h"
#include "h2lib/cluster.h"
#include "h2lib/hmatrix.h"
#include "FastBEM.h"

#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

#define DEBUG_ZMAT 0
#define FASTMETHOD 1

using namespace std;

int main(int argc, char * argv[])
{
    string filenameNode("./input/feko/total_nodes.txt");
    string filenameFace("./input/feko/all_surface.txt");
    string filenameEdge("./input/feko/surface_edge.txt");

	Mesh g_mesh(filenameNode, filenameEdge, filenameFace);
	g_mesh.display_mesh_info();
    BasisFunc g_basisfunc(g_mesh);

    double freq = 3e8;
    char integral_equation_type('p');

    #if !FASTMETHOD

    string filenameRCS("./output/b_rcs.csv");
    BasicBem bem(g_basisfunc);
    bem.build_euqation(freq, integral_equation_type);
    #if DEBUG_ZMAT
    bem.write_umat("output/b_exc.txt");
    bem.write_zmat("output/b_zmat.txt");
    #endif
    bem.solve();

    #else

    string filenameRCS("./output/f_rcs.csv");
    FastBemConf config = {100, 2, 1e-4};

	FastBem bem(g_basisfunc, config);

	auto t1 = Clock::now();
	bem.build_equation(freq, integral_equation_type);
	auto t2 = Clock::now();

	cout << "rel. error is :  " << scientific << bem.get_rel_error() << endl << endl;
	cout << "max rank is : " << bem.get_max_rank() << endl;
	cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0
         << " s" << endl;

//    OpenGLLib plot(bem.get_hmatrix());
//    plot.Display();

//    bem.print_zmat("./output/z.txt");
    cout << "print z over" << endl;

	bem.direct_solve_vec();

	bem.print_umat("./output/u.txt");
	cout << "print u over" << endl;
	auto t3 = Clock::now();
	cout << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() / 1000.0
         << " s" << endl;
	#endif // FASTMETHOD

	PostProcess post( g_basisfunc, bem.export_data(), bem.get_nrow(), bem.get_nrhs(), freq, integral_equation_type );
    post.gen_rcs(0, 360, 360);
    post.write_rcs(filenameRCS);

    return 0;
}
