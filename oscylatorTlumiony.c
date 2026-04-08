#include <stdio.h>
#include <math.h>

void sila(float* F, float *r, float *v, float k, float b, int n) 
{
    for(int i = 0; i < n; i++) 
	{
        F[i] = -k * r[i] - b * v[i];
    }
}

void euler(float *r, float *v, float m, float k, float b, float dt, int n) 
{
    float F[2];
    sila(F, r, v, k, b, n);
    for(int i = 0; i < n; i++) 
	{
        r[i] = r[i] + v[i] * dt;
        v[i] = v[i] + (F[i] / m) * dt;
    }
}

void verlet(float *r, float *v, float m, float k, float b, float dt, int n) 
{
    float Fk[2], Fk1[2], rk1[2], vk1_pred[2];
    int i;

    sila(Fk, r, v, k, b, n);

    for(i = 0; i < n; i++) 
	{
        rk1[i] = r[i] + v[i] * dt + (0.5 / m) * Fk[i] * dt * dt;
        vk1_pred[i] = v[i] + (Fk[i] / m) * dt;
    }

    sila(Fk1, rk1, vk1_pred, k, b, n);

    for(i = 0; i < n; i++) 
	{
        v[i] = v[i] + (0.5 / m) * (Fk[i] + Fk1[i]) * dt;
        r[i] = rk1[i];
    }
}

void symulacja(const char* nazwa, int tryb, float m, float k, float b, float dt, float t_max) 
{
    float r[2] = {0.5, 0.75}, r_old[2];;
    float v[2] = {0.0, 0.0};
    float t, ec, ec0, wop = 0;
    FILE *f = fopen(nazwa, "w");
    if(!f) return;
    
    ec0 = 0.5 * m * (v[0]*v[0] + v[1]*v[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);

    for(t = 0; t <= t_max; t += dt) 
	{
        ec = 0.5 * m * (v[0]*v[0] + v[1]*v[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);
        fprintf(f, "%f %f %f %f %f %f %f %f\n", t, r[0], r[1], v[0], v[1], ec, wop, ec- ec0);
        
        r_old[0] = r[0];
		r_old[1] = r[1];

        if(tryb == 0) euler(r, v, m, k, b, dt, 2);
        else verlet(r, v, m, k, b, dt, 2);
        
        wop += -b * (v[0] * (r[0] - r_old[0]) + v[1] * (r[1] - r_old[1]));
    }
    fclose(f);
}

int main() 
{
    float m = 0.5, k = 2.0, b = 0.1, dt = 0.01, t_max = 10.0;
    
    symulacja("euler.txt", 0, m, k, b, dt, t_max);
    symulacja("verlet.txt", 1, m, k, b, dt, t_max);
    
    return 0;
}
