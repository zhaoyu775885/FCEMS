#include "matrixplot.h"
#include <cmath>

pchmatrix g_hmatrix;

using namespace std;

const double half_gap = 1;

void OpenGLLib::InitGL(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

    glutCreateWindow("tst the opengl lib");
    glClearColor(0, 0, 0, 0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void
draw_red_block(pchmatrix hm, const double x, const double dx, const double y, const double dy)
{
    const double newx = 1100 -x;
    glBegin(GL_POLYGON);
    glColor3d(1,0,0);
    glVertex3f(newx-half_gap, y+half_gap, 0);
    glVertex3f(newx-dx+half_gap, y+half_gap, 0);
    glVertex3f(newx-dx+half_gap, y+dy-half_gap, 0);
    glVertex3f(newx-half_gap, y+dy-half_gap, 0);
    glEnd();
}

void
draw_green_block(pchmatrix hm, const double x, const double dx, const double y, const double dy)
{
    const double newx = 1100 -x;
    glBegin(GL_POLYGON);
    glColor3d(0,1,0);
    glVertex3f(newx-half_gap, y+half_gap, 0);
    glVertex3f(newx-dx+half_gap, y+half_gap, 0);
    glVertex3f(newx-dx+half_gap, y+dy-half_gap, 0);
    glVertex3f(newx-half_gap, y+dy-half_gap, 0);
    glEnd();
}

void draw_hmat(pchmatrix hm, const double beg_x, const double len_x, const double beg_y, const double len_y)
{
    if (hm->son == 0) {
        if (hm->f != 0) {
//            cout << "red ploted "
//                << beg_x << ", " << len_x << ", "
//                << beg_y << ", " << len_y << endl;
            draw_red_block(hm, beg_x, len_x, beg_y, len_y);
        }
        else if (hm->r != 0) {
            draw_green_block(hm, beg_x, len_x, beg_y, len_y);
        }
    }
    else {
        int nrow = hm->rc->size;
        int ncol = hm->cc->size;

        int prenrow = 0;
        for (int i=0;i<hm->rsons;++i) {
            int prencol = 0;
            for (int j=0;j<hm->csons;++j) {
                phmatrix newhm = hm->son[i+hm->rsons*j];
                double newx = beg_x + (float)prencol/ncol * len_x;
                double newlenx = (float)newhm->cc->size/ncol * len_x;
                double newy = beg_y + (float)prenrow/nrow * len_y;
                double newleny = (float)newhm->rc->size/nrow * len_y;
                draw_hmat(newhm, newx, newlenx, newy, newleny);
                prencol += newhm->cc->size;
            }
            prenrow += hm->son[i]->rc->size;
        }
    }

#if 0
    static double leaf = 100;
    if (len_x < leaf || len_y < leaf) {
        if (fabs(beg_x-beg_y) > 200) {
            draw_green_block(beg_x, len_x, beg_y, len_y);
//            cout << "draw green block" << endl;
        }
        else {
            draw_red_block(beg_x, len_x, beg_y, len_y);
//            cout << "draw red block" << endl;
        }
    }
    else {
        const double sub_len_x = len_x / 2.0;
        const double sub_len_y = len_y / 2.0;
        draw_hmat(beg_x, sub_len_x, beg_y, sub_len_y);
        draw_hmat(beg_x+sub_len_x, sub_len_x, beg_y, sub_len_y);
        draw_hmat(beg_x, sub_len_x, beg_y+sub_len_y, sub_len_y);
        draw_hmat(beg_x+sub_len_x, sub_len_x, beg_y+sub_len_y, sub_len_y);

    }
#endif // 0
}

void OpenGLLib::BlockDisplay( )
{
    glClear(GL_COLOR_BUFFER_BIT);
    glClear(GL_COLOR_BUFFER_BIT);

    const double origin_x = 50;
    const double origin_y = 50;

    const double canvas_x_len = 1000;
    const double canvas_y_len = 1000;

    glBegin(GL_POLYGON);
    glColor3d(1,1,1);
    glVertex3f(origin_x, origin_y, 0);
    glVertex3f(origin_x+canvas_x_len, origin_y, 0);
    glVertex3f(origin_x+canvas_x_len, origin_y+canvas_y_len, 0);
    glVertex3f(origin_x, origin_y+canvas_y_len, 0);
    glEnd();

    draw_hmat(g_hmatrix, origin_x, canvas_x_len, origin_y, canvas_y_len);

    /*
    glBegin(GL_POLYGON);
    glColor3d(0,1,0);
    glVertex3f(50.0, 0.0, 0);
    glVertex3f(50.0, 30.0, 0);
    glVertex3f(80.0, 30.0, 0);
    glVertex3f(80.0, 0.0, 0);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3d(1,1,1);
    glVertex3f(0.0, 10.0, 0);
    glVertex3f(10.0, 10.0, 0);
    glVertex3f(10.0, 0.0, 0);
    glVertex3f(0.0, 0.0, 0);
    glEnd();
    */

    glutSwapBuffers();
}

void OpenGLLib::HandleKeyPress(unsigned char key, int x, int y)
{
    cout << key << endl;
}

OpenGLLib::OpenGLLib(pchmatrix mat) : hmat(mat),
    win_xstart(0), win_ystart(0), win_width(WINDOW_WIDTH), win_height(WINDOW_HEIGHT)
{
    InitGL(0, 0);

    const double half_extend_width = EXTEND_WIDTH/2.0;
    const double half_extend_height = EXTEND_HEIGHT/2.0;

    gluOrtho2D(-half_extend_width, (GLdouble)win_width+half_extend_width, -half_extend_height, (GLdouble)win_height+half_extend_height);
}
void OpenGLLib::Display() const
{
    g_hmatrix = hmat;
    glutDisplayFunc(BlockDisplay);
    cout << "use opengl" << endl;
    glutKeyboardFunc(HandleKeyPress);
    glutMainLoop();
}

