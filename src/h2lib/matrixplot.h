#ifndef MATRIXPLOT_H
#define MATRIXPLOT_H

#include <GL/glut.h>
#include <cmath>
#include "hmatrix.h"

const double WINDOW_WIDTH(1920);
const double WINDOW_HEIGHT(1080);
const double EXTEND_WIDTH(0);
const double EXTEND_HEIGHT(0);

void
draw_green_block(pchmatrix hm, const double x, const double dx, const double y, const double dy);
void
draw_red_block(pchmatrix hm, const double x, const double dx, const double y, const double dy);
void draw_hmat(const double x, const double dx, const double y, const double dy);

class OpenGLLib
{
	public:
		OpenGLLib(pchmatrix mat);
        void Display() const;
        void InitGL(int argc, char **argv);
        static void BlockDisplay();
        static void HandleKeyPress(unsigned char key, int x, int y);


	private:
        pchmatrix hmat;
		double win_xstart;
		double win_ystart;
		double win_width;
		double win_height;
};

#endif // MATRIXPLOT_H
