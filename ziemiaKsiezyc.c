#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define G 6.67408e-11

double odleglosc(double *r1, double *r2, int n)
{
	int i;
	double r12 = 0;
	for(i = 0; i < n; i++)
	{	
		r12 += pow(r1[i] - r2[i], 2);
	}
	return sqrt(r12);
}
int sila(double *F, double M1, double *r1, double M2, double *r2, int n)
{
	int i;
	double Ftmp = -G * M1 * M2 / pow( odleglosc(r1,r2,n) ,3);
	for(i = 0; i < n; i++)
	{
		F[i] = Ftmp * (r1[i] - r2[i]);
	}
	return 0;
}

void euler(double *v1, double *v2, double M1, double *r1, double M2, double *r2, double dt, int n)
{
	double F1[2];
	double F2[2];
	
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

void verlet_predkosc(double *v1, double *v2, double M1, double *r1, double M2, double *r2, double dt, int n)
{
	double F1[2], F2[2];
	double F1_new[2], F2_new[2];
	double r1_new[2], r2_new[2];

	sila(F1, M1, r1, M2, r2, n);
	sila(F2, M2, r2, M1, r1, n);

	for(int i = 0; i < n; i++)
	{
		r1_new[i] = r1[i] + v1[i] * dt + 0.5 * (F1[i] / M1) * dt * dt;
		r2_new[i] = r2[i] + v2[i] * dt + 0.5 * (F2[i] / M2) * dt * dt;
	}

	sila(F1_new, M1, r1_new, M2, r2_new, n);
	sila(F2_new, M2, r2_new, M1, r1_new, n);

	for(int i = 0; i < n; i++)
	{
		v1[i] = v1[i] + 0.5 * ((F1[i] / M1) + (F1_new[i] / M1)) * dt;
		v2[i] = v2[i] + 0.5 * ((F2[i] / M2) + (F2_new[i] / M2)) * dt;

		r1[i] = r1_new[i];
		r2[i] = r2_new[i];
	}
}
void verlet_podst(double *v1, double *v2, double *r1, double *r2, double *r1_old, double *r2_old, double M1, double M2, double dt, int n)
{
	double F1[2], F2[2];
	double r1_new[2], r2_new[2];

	sila(F1, M1, r1, M2, r2, n);
	sila(F2, M2, r2, M1, r1, n);

	for(int i = 0; i < n; i++)
	{
		r1_new[i] = 2 * r1[i] - r1_old[i] + (F1[i] / M1) * dt * dt;
		r2_new[i] = 2 * r2[i] - r2_old[i] + (F2[i] / M2) * dt * dt;
		
		v1[i] = (r1_new[i] - r1_old[i]) / (2.0 * dt);
        v2[i] = (r2_new[i] - r2_old[i]) / (2.0 * dt);
	}

	for(int i = 0; i < n; i++)
	{
		r1_old[i] = r1[i];
		r2_old[i] = r2[i];

		r1[i] = r1_new[i];
		r2[i] = r2_new[i];
	}
}

void symulacja(const char* nazwa, int tryb, double mz, double mk, double dt, double t_max)
{
	double rz[2] = {0, 0}, rz_old[2];
	double rk[2] = {405696000, 0}, rk_old[2];
	double vz[2] = {0.0, 0.0};
	double vk[2] = {0, 968.0};
	double t, ec, ec0;
	double r_min = 1e12, r_max = 0;
	FILE *f = fopen(nazwa, "w");
	if(!f) return;

	ec0 = 0.5*mz*(vz[0]*vz[0]+vz[1]*vz[1]) + 0.5*mk*(vk[0]*vk[0]+vk[1]*vk[1]) - G*mz*mk/odleglosc(rz,rk,2);

	if(tryb == 2)
	{
		for(int i = 0; i < 2; i++) 
		{
		    rz_old[i] = rz[i];
		    rk_old[i] = rk[i];
		}
		
		double F1[2], F2[2];
		sila(F1, mz, rz, mk, rk, 2);
		sila(F2, mk, rk, mz, rz, 2);
		
		for(int i = 0; i < 2; i++) 
		{
		    rz[i] = rz[i] + vz[i] * dt + 0.5 * (F1[i] / mz) * dt * dt;
		    rk[i] = rk[i] + vk[i] * dt + 0.5 * (F2[i] / mk) * dt * dt;
		}
	}
	for(t = 0; t <= t_max; t += dt)
	{
		double d = odleglosc(rz, rk, 2);
        if(d < r_min) r_min = d;
        if(d > r_max) r_max = d;
        
		ec = 0.5*mz*(vz[0]*vz[0]+vz[1]*vz[1]) + 0.5*mk*(vk[0]*vk[0]+vk[1]*vk[1]) - G*mz*mk/odleglosc(rz,rk,2);
		
		double Px = mz * vz[0] + mk * vk[0];
        double Py = mz * vz[1] + mk * vk[1];
        double P_modul = sqrt(Px*Px + Py*Py);
        
        double x_rel = rk[0] - rz[0];
        double y_rel = rk[1] - rz[1];
        double vx_rel = vk[0] - vz[0];
        double vy_rel = vk[1] - vz[1];
		double iloczyn_wektorowy = fabs(x_rel * vy_rel - y_rel * vx_rel);
		
		double predkosc_polowa = 0.5 * iloczyn_wektorowy;
		
		fprintf(f, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n", t, rz[0], rz[1], rk[0], rk[1], vz[0], vz[1], vk[0], vk[1], ec, ec - ec0, P_modul, predkosc_polowa);

		if(tryb == 0) euler(vz, vk, mz, rz, mk, rk, dt, 2);
		else if(tryb == 1) verlet_predkosc(vz, vk, mz, rz, mk, rk, dt, 2);
		else verlet_podst(vz, vk, rz, rk, rz_old, rk_old, mz, mk, dt, 2);
	}
	fclose(f);
	
    double a = (r_min + r_max) / 2.0;
    double c = a - r_min;
    double b = sqrt(a * a - c * c);
    double e = c / a;

    printf("Wyniki dla: %s\n", nazwa);
    printf("a: %e, b: %e, c: %e, e: %f\n", a, b, c, e);
    printf("r_min (perygeum): %e, r_max (apogeum): %e\n\n", r_min, r_max);
}

int main()
{
	double mz = 5.97219e24, mk = 7.347673e22, dt = 360/2, t_max = 2592000;

	symulacja("euler.txt", 0, mz, mk, dt, t_max);
	symulacja("verlet_predkosc.txt", 1, mz, mk, dt, t_max);
	symulacja("verlet_podst.txt", 2, mz, mk, dt, t_max);

	return 0;
}// t polozenie ziemi polozenie ksiezyca predkosc ziemi predkosc ksiezyca
//  plot "euler.txt" u ($4-$2):($5-$3) w l
