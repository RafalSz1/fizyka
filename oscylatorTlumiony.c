#include <stdio.h>
#include <math.h>

int sila(float* F, float *r, float *v, float k, float b, int n) 
{
    for(int i = 0; i < n; i++) 
	{
        F[i] = -k * r[i] - b * v[i];
    }
    return 0;
}

int main() 
{
    float m = 0.5;
    float k = 2.0;
    float b = 0.1;
    float deltaT = 0.01;
    float t_max = 10.0;

    float r[2] = {0.5, 0.75};
    float v[2] = {0.0, 0.0};
    float F[2];

    float t, ek, ep, ec;
    FILE *plik = fopen("dane.txt", "w");

    if (plik == NULL) return 1;

	printf("t	 x 	  y	   Vx  	    Vy	     Ec\n\n");
	fprintf(plik, "t	 x 	  y	   Vx  	    Vy	     Ec\n\n");

    for (t = 0; t <= t_max; t += deltaT) 
	{

        ek = 0.5 * m * (v[0]*v[0] + v[1]*v[1]);
        ep = 0.5 * k * (r[0]*r[0] + r[1]*r[1]);
        ec = ek + ep;

        printf("%f\t%f\t%f\t%f\t%f\t%f\n", t, r[0], r[1], v[0], v[1], ec);
        fprintf(plik, "%f\t%f\t%f\t%f\t%f\t%f\n", t, r[0], r[1], v[0], v[1], ec);

        sila(F, r, v, k, b, 2);

        r[0] = r[0] + v[0] * deltaT;
        r[1] = r[1] + v[1] * deltaT;
        v[0] = v[0] + (F[0] / m) * deltaT;
        v[1] = v[1] + (F[1] / m) * deltaT;
    }

    fclose(plik);
    return 0;
}
