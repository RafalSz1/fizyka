#include <stdio.h>
#include <math.h>

#define G 6.67408e-11
#define MZ 5.972e24
#define MR 10000.0
#define RZ 6371000.0
#define RMAX 405696000.0
#define EPS_R 10.0
#define EPS_V 0.01

void symulacja(double alpha_deg, const char* filename) {
    FILE *f = fopen(filename, "w");
    double alpha = alpha_deg * M_PI / 180.0;
    double vII = sqrt(2.0 * G * MZ / RZ);

    double rx = RZ * cos(alpha), ry = RZ * sin(alpha);
    double vx = vII * cos(alpha), vy = vII * sin(alpha);
    double zx = 0, zy = 0, vzx = 0, vzy = 0;
    double t = 0, dt = 1.0;
    
    double Fold_x = 0, Fold_y = 0;

    while (sqrt(rx*rx + ry*ry) < RMAX) {
        double dx = rx - zx;
        double dy = ry - zy;
        double dist = sqrt(dx*dx + dy*dy);
        double dist3 = dist * dist * dist;

        double Fx = -G * MZ * MR * dx / dist3;
        double Fy = -G * MZ * MR * dy / dist3;

        double dtr = sqrt((2.0 * MR * EPS_R) / sqrt(Fx*Fx + Fy*Fy));
        double dtv = 1e6;
        
        if (t > 0) {
            double dFdt = sqrt(pow(Fx - Fold_x, 2) + pow(Fy - Fold_y, 2)) / dt;
            if (dFdt > 0) dtv = sqrt((2.0 * MR * EPS_V) / dFdt);
        }
        
        dt = (dtr < dtv) ? dtr : dtv;
        if (dt < 0.001) dt = 0.001;

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

int main() {
    symulacja(0.0, "wyniki_0.txt");
    symulacja(45.0, "wyniki_45.txt");
    return 0;
}
