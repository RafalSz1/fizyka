#include <stdio.h>
#include <math.h>

#define G 6.67408e-11
#define MZ 5.972e24
#define MK 7.342e22
#define MR 10000.0
#define RZ 6371000.0
#define RK 1737100.0
#define RMAX 405696000.0
#define VMIN 968.0
#define EPS_R 10.0
#define EPS_V 0.1

void symulacja(double alpha_deg, const char* filename) 
{
    FILE *f = fopen(filename, "w");
    double alpha = alpha_deg * M_PI / 180.0;
    double vII = sqrt(2.0 * G * MZ / RZ);

    double rZx = 0, rZy = 0;
    double vZx = 0, vZy = 0;
    
    double rKx = RMAX, rKy = 0;
    double vKx = 0, vKy = VMIN;
    double rRx = RZ * cos(alpha), rRy = RZ * sin(alpha);
    double vRx = vII * cos(alpha), vRy = vII * sin(alpha);
    double t = 0, dt = 1.0, dt_prev = 1.0;
    
    double Fold_x = 0, Fold_y = 0;
    
    printf("Symulacja dla alpha = %.2f stopni\n", alpha_deg);
    fprintf(f, "# t[s] rx[m] ry[m] kx[m] ky[m] E_total[J] dt[s] dtr dtv\n");

    while (1) 
	{
		double rRZ_x = rRx - rZx, rRZ_y = rRy - rZy;
        double rRK_x = rRx - rKx, rRK_y = rRy - rKy;
        double rKZ_x = rKx - rZx, rKZ_y = rKy - rZy;
		
		double dRZ = sqrt(rRZ_x * rRZ_x + rRZ_y * rRZ_y);
        double dRK = sqrt(rRK_x * rRK_x + rRK_y * rRK_y);
        double dKZ = sqrt(rKZ_x * rKZ_x + rKZ_y * rKZ_y);
        
        if (dRK <= RK) 
		{
            printf("SUKCES: Rakieta trafila w Ksiezyc! (t = %.2f h)\n", t / 3600.0);
            break;
        }
        if (dRK > RK && dRZ > dKZ) 
		{
            printf("PORAZKA: Rakieta minela Ksiezyc.\n");
            break;
        }

        double FRZx = -G * MR * MZ * rRZ_x / pow(dRZ, 3);
        double FRZy = -G * MR * MZ * rRZ_y / pow(dRZ, 3);
        double FRKx = -G * MR * MK * rRK_x / pow(dRK, 3);
        double FRKy = -G * MR * MK * rRK_y / pow(dRK, 3);

        double FRx = FRZx + FRKx;
        double FRy = FRZy + FRKy;
        double F_mag = sqrt(FRx * FRx + FRy * FRy);
        
        double dtr = sqrt((2.0 * MR * EPS_R) / F_mag);
        double dtv = 1e6;
        
        if (t > 0) 
		{			
			double dFx = FRx - Fold_x;
			double dFy = FRy - Fold_y;
			    
			double dF_dt_mag = sqrt(dFx*dFx + dFy*dFy) / dt;
			    
			if (dF_dt_mag > 0) 
			{
			    dtv = sqrt((2.0 * MR * EPS_V) / dF_dt_mag);
			}
        }
        
        dt_prev = dt;
        
        if (dtr < dtv) dt = dtr;
        else dt = dtv;
        
        if (dt > 1000.0) dt = 1000.0;
        if (dt < 0.01) dt = 0.01;

        Fold_x = FRx;
		Fold_y = FRy;
		
		double Ek = 0.5 * MZ * (vZx * vZx + vZy * vZy) +
                    0.5 * MK * (vKx * vKx + vKy * vKy) +
                    0.5 * MR * (vRx * vRx + vRy * vRy);
        double Ep = -G * (MZ * MK / dKZ + MZ * MR / dRZ + MK * MR / dRK);
        double E_total = Ek + Ep;
        
    	fprintf(f, "%e %e %e %e %e %e %e %e %e\n", t, rRx, rRy, rKx, rKy, E_total, dt, dtr, dtv);

		rRx += vRx * dt;
		rRy += vRy * dt;
        rKx += vKx * dt;
		rKy += vKy * dt;
        rZx += vZx * dt;
		rZy += vZy * dt;

		vRx += (FRx / MR) * dt;
        vRy += (FRy / MR) * dt;

		double FKZx = -G * MK * (MZ * rKZ_x / pow(dKZ, 3));
        double FKZy = -G * MK * (MZ * rKZ_y / pow(dKZ, 3));
        
        vKx += (FKZx - FRKx) / MK * dt;
        vKy += (FKZy - FRKy) / MK * dt;
        
        vZx += (-FKZx - FRZx) / MZ * dt;
        vZy += (-FKZy - FRZy) / MZ * dt;

        t += dt;
    }
    fclose(f);
}

int main() 
{
    symulacja(27, "trajektoria.txt");
    
    return 0;
}
