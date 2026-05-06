#include <stdio.h>
#include <math.h>

#define G 6.67408e-11
#define MZ 5.972e24
#define MR 10000.0
#define RZ 6371000.0
#define RMAX 405696000.0
#define EPS_R 10.0
#define EPS_V 0.1

void symulacja(double alpha_deg, const char* filename) 
{
    FILE *f = fopen(filename, "w");
    double alpha = alpha_deg * M_PI / 180.0;
    double vII = sqrt(2.0 * G * MZ / RZ);

    double rx = RZ * cos(alpha), ry = RZ * sin(alpha);
    double vx = vII * cos(alpha), vy = vII * sin(alpha);
    double zx = 0, zy = 0, vzx = 0, vzy = 0;
    double t = 0, dt = 1.0, dt_prev = 1.0;
    
    double Fold_x = 0, Fold_y = 0;

    while (sqrt(rx*rx + ry*ry) < RMAX) 
	{
        double dx = rx - zx;
        double dy = ry - zy;
        double dist = sqrt(dx*dx + dy*dy);
        double dist3 = dist * dist * dist;

        double Fx = -G * MZ * MR * dx / dist3;
        double Fy = -G * MZ * MR * dy / dist3;

        double dtr = sqrt((2.0 * MR * EPS_R) / sqrt(Fx*Fx + Fy*Fy));
        double dtv = 1000.0;
        
        if (t > 0) 
		{
			double dFx = Fx - Fold_x;
			double dFy = Fy - Fold_y;
			    
			double dF_dt_mag = sqrt(dFx*dFx + dFy*dFy) / dt;
			    
			if (dF_dt_mag > 0) 
			{
			    dtv = sqrt((2.0 * MR * dt_prev * EPS_V) / dF_dt_mag);
			}
        }
        
        dt_prev = dt;
        
        if (dtr < dtv) dt = dtr;
        else dt = dtv;
        
        if (dt > 1000.0) dt = 1000.0;
        if (dt < 0.01) dt = 0.01;

        fprintf(f, "%e %e %e %e %e\n", t, rx, ry, dtr, dtv);

        Fold_x = Fx;
		Fold_y = Fy;

        rx += vx * dt;
        ry += vy * dt;
        vx += (Fx / MR) * dt;
        vy += (Fy / MR) * dt;

        zx += vzx * dt;
        zy += vzy * dt;
        vzx += (-Fx / MZ) * dt;
        vzy += (-Fy / MZ) * dt;

        t += dt;
    }
    fclose(f);
    printf("Kat %.1f: Czas dolotu = %.2f h\n", alpha_deg, t / 3600.0);
}

int main() 
{
    symulacja(0.0, "wyniki_0.txt");
    symulacja(45.0, "wyniki_45.txt");
    return 0;
}
