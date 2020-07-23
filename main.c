#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.141592654

int iters = 2000;
int angle = 40;
double y_axis = 0;
double x_axis = 0;
double z_axis = 10000;

//产生符合正态分布的随机数
double gaussrand(double u, double z) {
    static float V1, V2, S;
    static int phase = 0;
    float X;
    if (phase == 0) {
        do {
            float U1 = (float) rand() / RAND_MAX;
            float U2 = (float) rand() / RAND_MAX;
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    phase = 1 - phase;
    return X*z + u;
}
//均匀分布采样
double _uniform(double min, double max, long int *seed) {
    double t = 0;
    *seed = 2045 * (*seed) + 1;
    *seed = *seed - (*seed / 1048576) * 1048576;
    t = (*seed) / 1048576.0;
    t = min + (max - min) * t;
    return t;
}


int main() {
    int numGrp = 46;
    double Paxis[46][3];
    double Daxis[46][3];
    double DL[46];
    double dx0 = z_axis * tan(angle * PI/180);
    double distL0[46];
    double distL1[46];
    double minus[46][3];
    double key[46];
    double Dxyz[3] = {0, 0, 0};
    for (int i = 4; i <= numGrp; i=i+2) {
        //变量初始化
        double calu[3] = {0, 0, 0};
        for (int k = 0; k < i; k++) {
            DL[k] = 0;
            distL0[k] = 0;
            distL1[k] = 0;
            key[k] = 0;
            for (int o = 0; o < 3; o++) {
                Paxis[k][o] = 0;
                Daxis[k][o] = 0;
                minus[k][o] = 0;
            }
        }
        for (int j = 1; j < iters; j++) {
            srand(time(0)); //确保每次随机不同
            double fun11 = 0, fun1 = 0, fun2 = 0, fun3 = 0, fun21 = 0, fun31 = 0; //fun1, 2, 3为方程，fun11，21，31为对应方程的d一阶导数
            double x = 0, y = 0, z = 0; //解
            
            int dX = 3;
            int dY = 3;
            int dZ = 5;
            int dL = 2;

            for (int k = 0; k < i; k++) {
                Daxis[k][0] = gaussrand(0, dX);
            }
            for (int k = 0; k < i; k++) {
                Daxis[k][1] = gaussrand(0, dY);
            }
            for (int k = 0; k < i; k++) {
                Daxis[k][2] = gaussrand(0, dZ);
            }
            for (int k = 0; k < i; k++) {
                DL[k] = gaussrand(0, dL);
            }
            for (int k = 0; k < i; k++) {
                Paxis[k][0] = rand()%2501 + dx0;
//                printf("%lf\n", Paxis[k][0]);
            }
            for (int k = 0; k < i; k++) {
                distL0[k] = sqrt(pow((x_axis - Paxis[k][0]), 2) + pow((y_axis - Paxis[k][1]), 2) + pow((z_axis - Paxis[k][2]), 2));
//                printf("%lf\n", distL0[k]);
            }
            for (int k = 0; k < i; k++) {
                distL1[k] = distL0[k] + DL[k];
            }
            for (int k = 0; k < i ; k++) {
                for (int p = 0; p < 3; p++) {
                    Paxis[k][p] += Daxis[k][p];
                }
            }
            //牛顿法求解非线性方程组
            for (int times = 0; times < 1000; times++) {
                for (int k = 0; k < i; k++) {
                    minus[k][0] = x_axis - Paxis[k][0];
                }
                for (int k = 0; k < i; k++) {
                    minus[k][1] = y_axis - Paxis[k][1];
                }
                for (int k = 0; k < i; k++) {
                    minus[k][2] = z_axis - Paxis[k][2];
                }
                for (int k = 0; k < i; k++) {
                    key[k] = pow(minus[k][0], 2) + pow(minus[k][1], 2) + pow(minus[k][2], 2) - pow(distL1[k], 2);
                }

                for (int k = 0; k < i; k++) {
                    fun1 += 4 * key[k] * minus[k][0];
                }
                for (int k = 0; k < i; k++) {
                    fun2 += 4 * key[k] * minus[k][1];
                }
                for (int k = 0; k < i; k++) {
                    fun3 += 4 * key[k] * minus[k][2];
                }
                for (int k = 0; k < i; k++) {
                    fun11 += 12 * pow(minus[k][0], 2) + 4 * pow(minus[k][1], 2) + 4 * pow(minus[k][2], 2);
                }
                for (int k = 0; k < i; k++) {
                    fun21 += 4 * pow(minus[k][0], 2) + 12 * pow(minus[k][1], 2) + 4 * pow(minus[k][2], 2);
                }
                for (int k = 0; k < i; k++) {
                    fun31 += 4 * pow(minus[k][0], 2) + 4 * pow(minus[k][1], 2) + 12 * pow(minus[k][2], 2);
                }
                x = x_axis - fun1/fun11;
                y = y_axis - fun2/fun21;
                z = z_axis - fun3/fun31;
                if (((x - x_axis) > 0?x - x_axis:x_axis - x) < 1e-5) {
                    break;
                }
                x_axis = x;
                y_axis = y;
                z_axis = z;
            }
            x_axis = 0;
            y_axis = 0;
            z_axis = 10000;

            Dxyz[0] = pow(x - x_axis, 2);
            Dxyz[1] = pow(y - y_axis, 2);
            Dxyz[2] = pow(z - z_axis, 2);
            for (int k = 0; k < 3; k++) {
                calu[k] += Dxyz[k];
            }
        }
        for (int k = 0; k < 3; k++) {
            calu[k] /= iters;
        }
        printf("num:%d   x:%lf; y:%lf; z:%lf\n", i, sqrt(calu[0]), sqrt(calu[1]), sqrt(calu[2]));
    }
    return 0;
}
