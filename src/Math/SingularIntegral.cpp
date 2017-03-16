#include "SingularIntegral.h"

void
Isp_and_Is(double *Isp, double *Is, const Point3D &pr1, const Point3D &pr2, const Point3D &pr3, const Point3D &pr0)
{
    double r1[3] = {pr1[0], pr1[1], pr1[2]};
    double r2[3] = {pr2[0], pr2[1], pr2[2]};
    double r3[3] = {pr3[0], pr3[1], pr3[2]};
    double r0[3] = {pr0[0], pr0[1], pr0[2]};

    double *vector_ri[3][2];
    double vector_n[3];
    double vector_ui[3][3];
    double scalar_Ri[3][2];
    double vector_Ii[3][3];
    double scalar_li[3][2];
    double scalar_d[3];
    double scalar_Pi0[3];
    double scalar_Ri02[3];
    double vector_temp[3];
    double vector_Pi0[3][3];

    vector_ri[0][0] = r1;
    vector_ri[0][1] = r2;
    vector_ri[1][0] = r2;
    vector_ri[1][1] = r3;
    vector_ri[2][0] = r3;
    vector_ri[2][1] = r1;

    Isp[0] = Isp[1] = Isp[2] = 0;
    *Is = 0;

    double ab[3], bc[3];
    vectorSub(vector_ri[0][1], vector_ri[0][0], ab, 3);
    vectorSub(vector_ri[1][1], vector_ri[1][0], bc, 3);
    crossProduct(ab, bc, vector_n);
    normalizeVector(vector_n, 3);
    for (int i=0;i<3;i++) {
        vectorSub(vector_ri[i][1], vector_ri[i][0], vector_Ii[i], 3);
        normalizeVector(vector_Ii[i], 3);
        crossProduct(vector_Ii[i], vector_n, vector_ui[i]);
        scalar_Ri[i][0] = lenEdge(vector_ri[i][0], r0, 3);
        scalar_Ri[i][1] = lenEdge(vector_ri[i][1], r0, 3);
        vectorSub(vector_ri[i][0], r0, vector_temp, 3);
        scalar_li[i][0] = dotProduct(vector_Ii[i], vector_temp, 3);
        vectorSub(vector_ri[i][1], r0, vector_temp, 3);
        scalar_li[i][1] = dotProduct(vector_Ii[i], vector_temp, 3);
        vectorSub(r0, vector_ri[i][0], vector_temp, 3);
        scalar_d[i] = dotProduct(vector_temp, vector_n, 3);
        scalar_Pi0[i] = fabs( dotProduct(vector_temp, vector_ui[i], 3) );
        scalar_Ri02[i] = scalar_Pi0[i]*scalar_Pi0[i] + scalar_d[i]*scalar_d[i];
        for (int j=0;j<3;j++) {
            vector_Pi0[i][j] = (vector_ri[i][1][j] - r0[j] - scalar_li[i][1] * vector_Ii[i][j])/scalar_Pi0[i];
        }
        for (int j=0;j<3;j++) {
            Isp[j] += (1.0/2)*(scalar_Ri02[i]*log((scalar_Ri[i][1]+scalar_li[i][1])/(scalar_Ri[i][0]+scalar_li[i][0])) + \
            scalar_li[i][1]*scalar_Ri[i][1] - scalar_li[i][0]*scalar_Ri[i][0]) * vector_ui[i][j];
        }
        *Is +=  dotProduct(vector_Pi0[i], vector_ui[i], 3) * (scalar_Pi0[i] * log((scalar_Ri[i][1]+scalar_li[i][1])/ \
        (scalar_Ri[i][0]+scalar_li[i][0])) - fabs(scalar_d[i])*(atan(scalar_Pi0[i]*scalar_li[i][1]/(scalar_Ri02[i]+\
        abs(scalar_d[i])*scalar_Ri[i][1])) - atan(scalar_Pi0[i]*scalar_li[i][0]/(scalar_Ri02[i]+abs(scalar_d[i])*
        scalar_Ri[i][0]))));
    }
}


void
Isp3_and_Isp(double *Isp3, double *Isp,
            const Point3D &pr1, const Point3D &pr2, const Point3D &pr3, const Point3D &pr0)
{
    double r1[3] = {pr1[0], pr1[1], pr1[2]};
    double r2[3] = {pr2[0], pr2[1], pr2[2]};
    double r3[3] = {pr3[0], pr3[1], pr3[2]};
    double r0[3] = {pr0[0], pr0[1], pr0[2]};
    double *vector_ri[3][2];
    double vector_n[3];
    double vector_ui[3][3];
    double scalar_Ri[3][2];
    double vector_Ii[3][3];
    double scalar_li[3][2];
    double scalar_d[3];
    double scalar_Pi0[3];
    double scalar_Ri02[3];
    double vector_temp[3];

    vector_ri[0][0] = r1;
    vector_ri[0][1] = r2;
    vector_ri[1][0] = r2;
    vector_ri[1][1] = r3;
    vector_ri[2][0] = r3;
    vector_ri[2][1] = r1;

    Isp3[0] = Isp3[1] = Isp3[2] = 0;
    Isp[0] = Isp[1] = Isp[2] = 0;

    double ab[3], bc[3];
    vectorSub(vector_ri[0][1], vector_ri[0][0], ab, 3);
    vectorSub(vector_ri[1][1], vector_ri[1][0], bc, 3);
    crossProduct(ab, bc, vector_n);
    normalizeVector(vector_n, 3);
    for (int i=0;i<3;i++) {
    vectorSub(vector_ri[i][1], vector_ri[i][0], vector_Ii[i], 3);
    normalizeVector(vector_Ii[i], 3);
    crossProduct(vector_Ii[i], vector_n, vector_ui[i]);
    scalar_Ri[i][0] = lenEdge(vector_ri[i][0], r0, 3);
    scalar_Ri[i][1] = lenEdge(vector_ri[i][1], r0, 3);
    vectorSub(vector_ri[i][0], r0, vector_temp, 3);
    scalar_li[i][0] = dotProduct(vector_Ii[i], vector_temp, 3);
    vectorSub(vector_ri[i][1], r0, vector_temp, 3);
    scalar_li[i][1] = dotProduct(vector_Ii[i], vector_temp, 3);
    vectorSub(r0, vector_ri[i][0], vector_temp, 3);
    scalar_d[i] = dotProduct(vector_temp, vector_n, 3);
    //printf("scalar_d: %e\n",scalar_d[i]);
    scalar_Pi0[i] = fabs( dotProduct(vector_temp, vector_ui[i], 3) );
    scalar_Ri02[i] = scalar_Pi0[i]*scalar_Pi0[i] + scalar_d[i]*scalar_d[i];

    for (int j=0;j<3;j++) {
        Isp[j] += (1.0/2)*(scalar_Ri02[i]*log((scalar_Ri[i][1]+scalar_li[i][1])/(scalar_Ri[i][0]+scalar_li[i][0])) + \
        scalar_li[i][1]*scalar_Ri[i][1] - scalar_li[i][0]*scalar_Ri[i][0]) * vector_ui[i][j];
        Isp3[j] -= log((scalar_Ri[i][1]+scalar_li[i][1])/(scalar_Ri[i][0]+scalar_li[i][0])) * vector_ui[i][j];
        }
    }
}



/*
void sinItg(halfrtpBF hrbf, double rm[3], double *singlar1, double *singlar2)
{
    double a, b, c, d;
    double av[3];

//    double tmp[3];
//    vectorSub(hrbf.v1, hrbf.v4, tmp, 3);
//    normalizeVector(tmp, 3);
    vectorSub(hrbf.v2, hrbf.v1, av, 3);
//    numProduct(av, dotProduct(tmp, hrbf.vct, 3), 3);
    normalizeVector(av, 3);

    double tmp1[3], tmp2[3], tmp3[3], tmp4[4];
    vectorSub(hrbf.v1, rm, tmp1, 3);
    vectorSub(hrbf.v2, rm, tmp2, 3);
    vectorSub(hrbf.v3, rm, tmp3, 3);
    vectorSub(hrbf.v4, rm, tmp4, 3);

    a = dotProduct(tmp1, hrbf.vct, 3);
    b = dotProduct(tmp4, hrbf.vct, 3);
    c = dotProduct(tmp3, av, 3);
    d = dotProduct(tmp4, av, 3);

    a = a<b?a:b;
    b = hrbf.lx+a;
    c = c<d?c:d;
    d = hrbf.ly + c;

    *singlar1 = sinItgFml1(b, d)+sinItgFml1(a, c)-sinItgFml1(b, c)-sinItgFml1(a, d);
    *singlar2 = sinItgFml2(b, d)+sinItgFml2(a, c)-sinItgFml2(b, c)-sinItgFml2(a, d);
}

double sinItgFml2(double x, double y)
{
    return ( x * log(y+sqrt(x*x+y*y)) + y * log(x+sqrt(x*x+y*y)) );
}

double sinItgFml1(double x, double y)
{
    return 0.5 * ( y*sqrt(x*x+y*y) + x*x*log(y+sqrt(x*x+y*y)) );
}
*/
