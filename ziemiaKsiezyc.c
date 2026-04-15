#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define G 6.67408e-11

float odleglosc(float *r1, float *r2, int n)
{
	int i;
	float r12 = 0;
	for(i = 0; i < n; i++)
	r12 += pow(r1[i] - r2[i], 2);
	return sqrt(r12);
}
int sila(float *F, float M1, float *r1, float M2, float *r2, int n)
{
	int i;
	float Ftmp = -G * M1 * M2 / pow( odleglosc(r1,r2,n) , 3 );
	for(i = 0; i < n; i++)
	F[i] = Ftmp * ( r1[i] - r2[i] );
	return 0;
}

void euler(float *v1, float *v2, float M1, float *r1, float M2, float *r2, float dt, int n)
{
	float F1[2];
	float F2[2];
	sila(F1, M1, r1, M2, r2, n);
	sila(F2, M2, r2, M1, r1, n);
	for(int i = 0; i < n; i++)
	{
		r1[i] = r1[i] + v1[i] * dt;
		v1[i] = v1[i] + (F1[i] / M1) * dt;
		
		r2[i] = r2[i] + v2[i] * dt;
		v2[i] = v2[i] + (F2[i] / M2) * dt;
	}
}

//void verlet_predkosc(float *r, float *v, float m, float k, float dt, int n)
//{
//	float Fk[2], Fk1[2], rk1[2], vk1_pred[2];
//	int i;
//
//	sila(Fk, r, k, n);
//
//	for(i = 0; i < n; i++)
//	{
//		rk1[i] = r[i] + v[i] * dt + (0.5 / m) * Fk[i] * dt * dt;
//	}
//
//	sila(Fk1, rk1, k, n);
//
//	for(i = 0; i < n; i++)
//	{
//		v[i] = v[i] + (0.5 / m) * (Fk[i] + Fk1[i]) * dt;
//		r[i] = rk1[i];
//	}
//}
//
//void verlet_podst(float *r,float *v, float *r_old, float m, float k, float dt, int n)
//{
//	float F[2], r_new[2];
//	int i;
//
//	sila(F, r, k, n);
//
//	for(i = 0; i < n; i++)
//	{
//		r_new[i] = 2*r[i] - r_old[i] + (F[i]/m)*dt*dt;
//	}
//	for(i = 0; i < n; i++)
//    {
//        v[i] = (r_new[i] - r_old[i]) / (2*dt);
//    }
//	for(i = 0; i < n; i++)
//	{
//		r_old[i] = r[i];
//		r[i] = r_new[i];
//	}
//}

void symulacja(const char* nazwa, int tryb, float mz, float mk, float dt, float t_max)
{
	float rz[2] = {0, 0}, rz_old[2];
	float rk[2] = {405696000, 0}, rk_old[2];
	float vz[2] = {0.0, 0.0};
	float vk[2] = {0, 968.0};
	float t, ec, ec0;
	FILE *f = fopen(nazwa, "w");
	if(!f) return;

	//ec0 = 0.5 * m * (v[0]*v[0] + v[1]*v[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);

//	if(tryb == 2)
//	{
//	    float F[2];
//	    sila(F, r, k, 2);
//	
//	    for(int i = 0; i < 2; i++)
//	    {
//	        r_old[i] = r[i] - v[i]*dt + (F[i]/(2*m))*dt*dt;
//	    }
//	}
	for(t = 0; t <= t_max; t += dt)
	{
		//ec = 0.5 * m * (v[0]*v[0] + v[1]*v[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);
		fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f\n", t, rz[0], rz[1], rk[0], rk[1], vz[0], vz[1], vk[0], vk[1], ec, ec - ec0);

		if(tryb == 0) euler(vz, vk, mz, rz, mk, rk, dt, 2);
		//else if(tryb == 1) verlet_predkosc(r, v, m, k, dt, 2);
		//else verlet_podst(r, v, r_old, m, k, dt, 2);
	}
	fclose(f);
}

int main()
{
	float mz = 5.97219e24, mk = 7.347673e22, dt = 360, t_max = 2592000;

	symulacja("euler.txt", 0, mz, mk, dt, t_max);
	symulacja("verlet_predkosc.txt", 1, mz, mk, dt, t_max);
	symulacja("verlet_podst.txt", 2, mz, mk, dt, t_max);

	return 0;
}// t polozenie ziemi polozenie ksiezyca predkosc ziemi predkosc ksiezyca
