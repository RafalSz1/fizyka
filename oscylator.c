#include <stdio.h>
#include <math.h>

void sila(float* F, float *r, float k, int n)
{
	for(int i = 0; i < n; i++)
	{
		F[i] = -k * r[i];
	}
}

void euler(float *r, float *v, float m, float k, float dt, int n)
{
	float F[2];
	sila(F, r, k, n);
	for(int i = 0; i < n; i++)
	{
		r[i] = r[i] + v[i] * dt;
		v[i] = v[i] + (F[i] / m) * dt;
	}
}

void verlet_predkosc(float *r, float *v, float m, float k, float dt, int n)
{
	float Fk[2], Fk1[2], rk1[2], vk1_pred[2];
	int i;

	sila(Fk, r, k, n);

	for(i = 0; i < n; i++)
	{
		rk1[i] = r[i] + v[i] * dt + (0.5 / m) * Fk[i] * dt * dt;
	}

	sila(Fk1, rk1, k, n);

	for(i = 0; i < n; i++)
	{
		v[i] = v[i] + (0.5 / m) * (Fk[i] + Fk1[i]) * dt;
		r[i] = rk1[i];
	}
}

void verlet_podst(float *r,float *v, float *r_old, float m, float k, float dt, int n)
{
	float F[2], r_new[2];
	int i;

	sila(F, r, k, n);

	for(i = 0; i < n; i++)
	{
		r_new[i] = 2*r[i] - r_old[i] + (F[i]/m)*dt*dt;
	}
    
	for(i = 0; i < n; i++)
	{
		r_old[i] = r[i];
		r[i] = r_new[i];
	}
}

void symulacja(const char* nazwa, int tryb, float m, float k, float dt, float t_max)
{
	float r[2] = {0.5, 0.75}, r_old[2];
	float v[2] = {0.0, 0.0};
	float t, ec, ec0;
	FILE *f = fopen(nazwa, "w");
	if(!f) return;

	ec0 = 0.5 * m * (v[0]*v[0] + v[1]*v[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);

	if(tryb == 2)
	{
	    float r_old[2], r_new[2], F[2], v_ec[2];
	    int i;
	    sila(F, r, k, 2);
	
	    for(int i = 0; i < 2; i++)
	    {
	        r_old[i] = r[i] - v[i]*dt + (F[i]/(2*m))*dt*dt;
	    }
	    for(i = 0; i < 2; i++)
        {
			r_new[i] = 2*r[i] - r_old[i] + (F[i]/m)*dt*dt;
		}
		for(t = 0; t <= t_max; t += dt)
        {
            for(i = 0; i < 2; i++)
            {
            	v_ec[i] = (r_new[i] - r_old[i]) / (2 * dt);
			}

            ec = 0.5 * m * (v_ec[0]*v_ec[0] + v_ec[1]*v_ec[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);

            fprintf(f, "%f %f %f %f %f %f %f\n", t, r[0], r[1], v_ec[0], v_ec[1], ec, ec - ec0);

            if(t + dt > t_max) break;
            
            for(i = 0; i < 2; i++)
            {
                r_old[i] = r[i];
                r[i] = r_new[i];
            }

            sila(F, r, k, 2);
            for(i = 0; i < 2; i++)
            {
            	r_new[i] = 2*r[i] - r_old[i] + (F[i]/m)*dt*dt;
			}
        }
        fclose(f);
	}
	for(t = 0; t <= t_max; t += dt)
	{
		ec = 0.5 * m * (v[0]*v[0] + v[1]*v[1]) + 0.5 * k * (r[0]*r[0] + r[1]*r[1]);
		fprintf(f, "%f %f %f %f %f %f %f\n", t, r[0], r[1], v[0], v[1], ec, ec - ec0);

		if(tryb == 0) euler(r, v, m, k, dt, 2);
		else if(tryb == 1) verlet_predkosc(r, v, m, k, dt, 2);
	}
	fclose(f);
}

int main()
{
	float m = 0.5, k = 2.0, dt = 0.01, t_max = 10.0;

	symulacja("euler.txt", 0, m, k, dt, t_max);
	symulacja("verlet_predkosc.txt", 1, m, k, dt, t_max);
	symulacja("verlet_podst.txt", 2, m, k, dt, t_max);

	return 0;
}
