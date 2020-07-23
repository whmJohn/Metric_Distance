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
double _uniform(double min, double max, long int seed) {
    double t = 0;
    seed = 2045 * (seed) + 1;
    seed = seed - (seed / 1048576) * 1048576;
    t = (seed) / 1048576.0;
    t = min + (max - min) * t;
    return t;
}

//牛顿法求非线性方程组
double * noneliear(int num, double x_axis, double y_axis, double z_axis, double Paxis[num][3], double distL1[num], double sol4_fsolve[num]) {
    double minus[num][3];
    double fun11 = 0, fun1 = 0, fun2 = 0, fun3 = 0, fun21 = 0, fun31 = 0; //fun1, 2, 3为方程，fun11，21，31为对应方程的d一阶导数
    double x = 0, y = 0, z = 0; //解
    double key[num];
    for (int times = 0; times < 1000; times++) {
        for (int k = 0; k < num; k++) {
            minus[k][0] = x_axis - Paxis[k][0];
            key[k] = 0;
        }
        for (int k = 0; k < num; k++) {
            minus[k][1] = y_axis - Paxis[k][1];
        }
        for (int k = 0; k < num; k++) {
            minus[k][2] = z_axis - Paxis[k][2];
        }
        for (int k = 0; k < num; k++) {
            key[k] = pow(minus[k][0], 2) + pow(minus[k][1], 2) + pow(minus[k][2], 2) - pow(distL1[k], 2);
        }

        for (int k = 0; k < num; k++) {
            fun1 += 4 * key[k] * minus[k][0];
        }
        for (int k = 0; k < num; k++) {
            fun2 += 4 * key[k] * minus[k][1];
        }
        for (int k = 0; k < num; k++) {
            fun3 += 4 * key[k] * minus[k][2];
        }
        for (int k = 0; k < num; k++) {
            fun11 += 12 * pow(minus[k][0], 2) + 4 * pow(minus[k][1], 2) + 4 * pow(minus[k][2], 2);
        }
        for (int k = 0; k < num; k++) {
            fun21 += 4 * pow(minus[k][0], 2) + 12 * pow(minus[k][1], 2) + 4 * pow(minus[k][2], 2);
        }
        for (int k = 0; k < num; k++) {
            fun31 += 4 * pow(minus[k][0], 2) + 4 * pow(minus[k][1], 2) + 12 * pow(minus[k][2], 2);
        }
        x = x_axis - fun1/fun11;
        y = y_axis - fun2/fun21;
        z = z_axis - fun3/fun31;
        if (((x - x_axis) > 0?x - x_axis:x_axis - x) < 1e-5) {
            sol4_fsolve[0] = x;
            sol4_fsolve[1] = y;
            sol4_fsolve[2] = z;
            break;
        }
        x_axis = x;
        y_axis = y;
        z_axis = z;
    }
    return sol4_fsolve;
}

int main() {
    //78
    int num = 4;
    double axis[num][3];
    double y_axis = _uniform(-20000, 20000, 102);
    // printf("yaxis%lf", y_axis);
    double z_axis = _uniform(6000, 18000, 102); 
    double x_axis = _uniform(-10000, 10000, 102); 
    // printf("xaxis:%lf", x_axis);
    double slope_distance[num];
    for (int i = 0; i < num; i++) {
        axis[i][1] = y_axis;
        axis[i][2] = _uniform(0, 4000, 102+10*i); 
        axis[i][0] = _uniform(6000, 8500, 102+10*i);  
    }
    double b[num];
    for (int i = 0; i < num; i++) {
        b[i] = fabs(z_axis - axis[i][2]);
        // printf("%lf", b[i]);
    } 
    for (int i = 0; i < num; i++) {
        slope_distance[i] = sqrt(pow(x_axis - axis[i][0], 2) + pow(z_axis - axis[i][2], 2));
        // printf("slopedis;%lf", slope_distance[i]);
    }
    //79
    for (int i = 0; i < num; i++) {
        printf("axisx:%lf, y:%lf, z:%lf\n", axis[i][0], axis[i][1], axis[i][2]);
    }
    //80
    int iter = 100;
    double calu[3];
    double range_xl[17];
    double gussian_noise[num];
    double dectected[num];
    double * solve;
    double sol4_fsolve[num];
    double results[17][3];
    for (int i = 0; i < 17; i++) {
        range_xl[i] = i * 0.5 + 1;
    }
    for (int i = 0; i < 17; i++) {
        for (int j = 0; j < iter; j++) {
            srand(time(0)); //确保每次随机不同
            for (int o = 0; o < num; o++) {
                gussian_noise[o] = gaussrand(0, range_xl[i]);
                dectected[o] = gussian_noise[o] + slope_distance[o];
            }

            solve = noneliear(num, x_axis, y_axis, z_axis, axis, dectected, sol4_fsolve);//Paxis以及斜距确定*****
            // printf("solve0:%lf,solve1:%lf,solve2:%lf\n", solve[0], solve[1], solve[2]);
            calu[0] += pow(solve[0] - x_axis, 2);
            calu[1] += pow(solve[1] - y_axis, 2);
            calu[2] += pow(solve[2] - z_axis, 2);
        }
        for (int q = 0; q < 3; q++) {
            calu[q] /= iter;
            results[i][q] = sqrt(calu[q]);
            printf("results:%lf,", results[i][q]);
            if(q == 2) {
                printf("\n");
            }
        }
    }
    //21
    int count = 100;
    double results1[count][10][3];
    double fresults[10][3];
    while(count) {
        int myint = rand();
        for (int k = 3; k < 13; k++) {
            num = k;
            srand(myint);
            double axis1[num][3];
            y_axis = _uniform(-20000, 20000, myint);
            z_axis = _uniform(6000, 18000, myint);
            x_axis = _uniform(-10000, 10000, myint);
            for (int i = 0; i < num; i++) {
                axis1[i][1] = y_axis;
                axis1[i][2] = _uniform(0, 4000, myint + 10*i); 
                axis1[i][0] = _uniform(x_axis+10000, x_axis+12500, myint + 10*i); 
            }
            double slope_distance1[num];
            for (int i = 0; i < num; i++) {
                slope_distance1[i] = sqrt(pow(x_axis - axis1[i][0], 2) + pow(z_axis - axis1[i][2], 2));
            }
            double sigma = 3.0;
            double gussian_noise1[num];
            double dectected1[num];
            double sol4_fsolve1[num];
            for (int i = 0; i < iter; i++) {
                srand(time(0)); //确保每次随机不同
                for (int o = 0; o < num; o++) {
                    gussian_noise1[o] = gaussrand(0, sigma);
                    dectected1[o] = gussian_noise1[o] + slope_distance1[o];
                }
                solve = noneliear(num, x_axis, y_axis, z_axis, axis1, dectected1, sol4_fsolve1);//Paxis以及斜距确定
                calu[0] += pow(solve[0] - x_axis, 2);
                calu[1] += pow(solve[1] - y_axis, 2);
                calu[2] += pow(solve[2] - z_axis, 2);
            }
            for (int q = 0; q < 3; q++) {
                calu[q] /= iter;
                results1[100-count][k-3][q] = sqrt(calu[q]);
            }
        }
        count -= 1;
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 3; j++) {
            for (int q = 0; q < 100; q++) {
                fresults[i][j] += results1[q][i][j]; 
            }
            fresults[i][j] /= 100;
            printf("final_result:%lf", fresults[i][j]);
        }
    }

    //9
    count = 100;
    double results2[count][4][3];
    double fresults2[4][3];
    while(count) {
        int myint = rand();
        num = 12;
        srand(myint);
        double axis2[num][3];
        y_axis = _uniform(-20000, 20000, myint);
        z_axis = _uniform(8000, 18000, myint);
        x_axis = _uniform(-10000, 10000, myint);
        for (int i = 0; i < num; i++) {
            axis2[i][1] = y_axis;
            axis2[i][2] = _uniform(0, 2000, myint + 10*i); 
            axis2[i][0] = _uniform(x_axis+10000, x_axis+12500, myint + 10*i); 
        }
        double gt[num][3];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                if (j == 0) {
                    gt[i][j] = x_axis;
                }
                if (j == 1) {
                    gt[i][j] = y_axis + 1.5*i; 
                }
                if (j == 2) {
                    gt[i][j] = z_axis;
                }
            }
        }
        double slope_distance2[num];
        for (int i = 0; i < num; i++) {
            slope_distance2[i] = sqrt(pow(x_axis - axis2[i][0], 2) + pow(z_axis - axis2[i][2], 2));
        }
        double sigma = 3.0;
        double gussian_noise2[num];
        double dectected2[num];
        double sol4_fsolve2[num];
        double calu2[4][3];
        for (int i = 0; i < iter; i++) {
            srand(time(0)); //确保每次随机不同
            for (int o = 0; o < num; o++) {
                gussian_noise2[o] = gaussrand(0, sigma);
                dectected2[o] = gussian_noise2[o] + slope_distance2[o];
            }
            solve = noneliear(num, x_axis, y_axis, z_axis, axis2, dectected2, sol4_fsolve2);//Paxis以及斜距确定
            for (int i = 0 ;i < 4; i++) {
                for (int j = 0; j < 3; j++) {
                    calu2[i][j] += pow(solve[0] - gt[i][j], 2);
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            for (int q = 0; q < 3; q++) {
                calu2[i][q] /= iter;
                results2[100-count][i][q] = sqrt(calu2[i][q]);
            }
        } 
        count -= 1;
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 3; j++) {
            for (int q = 0; q < 100; q++) {
                fresults2[i][j] += results2[q][i][j]; 
            }
            fresults2[i][j] /= 100;
            printf("final_result2:%lf", fresults[i][j]);
        }
    }

    //31 39 35 同80
        
    


    return 0;
}
