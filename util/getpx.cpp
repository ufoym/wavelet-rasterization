#include <math.h>
#define DLLEXPORT extern "C" __declspec(dllexport)
DLLEXPORT void get_px(float area, int max_j, int px_num, 
                     float * c_arr, float * lattice_arr, float * px_arr)
{
    float Ex[] = {0, 1, 1};
    float Ey[] = {1, 0, 1};
    float * p_lattice = lattice_arr;
    for (int i = 0; i < px_num; ++i) {
        float px = *p_lattice++;
        float py = *p_lattice++;
        float s = area;
        float * pc = c_arr;
        for (int j = 0; j <= max_j; ++j) {
            float exp_j = pow(2, j);
            float exp_jpx = exp_j * px;
            float exp_jpy = exp_j * py;
            for (int kx = 0; kx < exp_j; ++kx) {
                for (int ky = 0; ky < exp_j; ++ky) {
                    float exp_jpkx = exp_jpx-kx;
                    float exp_jpky = exp_jpy-ky;
                    for (int i = 0; i < 3; ++i) {
                        float c = *pc++;
                        if (exp_jpkx < 0 || exp_jpkx >= 1)  continue;
                        if (exp_jpky < 0 || exp_jpky >= 1)  continue;
                        bool neg_x = (exp_jpkx >= 0.5 && Ex[i] != 0);
                        bool neg_y = (exp_jpky >= 0.5 && Ey[i] != 0);
                        if (neg_x && (!neg_y) || (!neg_x) && neg_y)
                            s -= c;
                        else  
                            s += c;
                    }
                }
            }
        }
        px_arr[i] = s;       
    }
}