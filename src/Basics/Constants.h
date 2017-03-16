#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#include <complex>

namespace Constants
{
    const double PI = 3.1415926535897932384626433832795029;
    const double M_PI1 = PI;
    const double PI2 = 6.2831853071795864769252867665590058;
    const double M_PI_21 = 1.5707963267948966192313216916397514;/* pi/2*/
    const double M_2_PI1 = 0.63661977236758134308;	/* 2/pi */

    const double SQRT2 = 1.41421356237309504880;
    const double SQRT3 = 1.73205080756887729353;

    const double EPSILON0 = 8.854187817e-12;
    const double M_epsilon = 8.8541878176e-15; // (unit: F/mm)
    const double MU0 = 4.0e-7*PI;
    const double M_mu = M_PI1*4.0e-10; // (unit: H/mm)
    const double SPEED_OF_LIGHT = 2.99792458e8;
    const double ETA0 = 376.7303135;
    const double M_Z0 = M_PI1*119.9169832; // free space impedance
    const double M_GAMMA = 0.5772;

    const double INF = 1e16;
    const double TOLERANCE = 1e-4;
    const double SOLVER_TOLERANCE = 1e-12;

    const std::complex<double> CJ(0.,1.0);

//    extern const char* MatlabColor;
//
//    //some defaults
//    extern int NUM_SEGMENTS_CIRCLE;
//
//    // product related
//    const std::string STACKUP_EXTENSION=".lyr";
//    const std::string PARAMAX_EXTENSION=".pmax";
//    const std::string RLCGDATA_EXTENSION=".rlcg";
//
//    const std::string TLINE2D_PRODUCT_NAME ="TLine2D";
//	const std::string MESH_EXTENSION=".mesh";
//    const std::string CAD_EXTENSION=".geom";
//	const std::string QUICKPOWER_EXTENSION = ".esf";
//
//	const double Cpmax   = 1.0e+6; // Threshold value of decoupling capacitance
//
//	const double xk[]={-0.861136311594, -0.339981043585, 0.339981043585,  0.861136311594};
//	const double wk[]={0.347854845137, 0.652145154863, 0.652145154863, 0.347854845137};
/*
	const double xk[]={-0.9894009349, -0.9445750230,\
		-0.8656312023, -0.7554044083,\
		-0.6178762444, -0.4580167776,\
		-0.2816035507, -0.09501250983,\
		0.09501250983,  0.2816035507,\
		0.4580167776,  0.6178762444,\
		0.7554044083,  0.8656312023,\
		0.9445750230,  0.9894009349};

	const double wk[]={0.02715245941, 0.06225352393,\
		0.09515851168,  0.1246289712,\
		0.1495959888,  0.1691565193,\
		0.1826034150,  0.1894506104,\
		0.1894506104,  0.1826034150,\
		0.1691565193,  0.1495959888,\
		0.1246289712, 0.09515851168,\
		0.06225352393, 0.02715245941};
*/
	const double g_dDefaultLineLengthForRLGC = 1e-4;
};

#endif
