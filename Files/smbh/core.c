#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

#define CONCENTRATION_PARAMETER 4
#define MATTER_DENSITY_PARAMETER 0

volatile double R_VIR = 0;
volatile double HALO_MASS = 1e3;
volatile double SMBH_MASS = 1;
volatile double BARYONIC_TOTAL_MASS = 158; //1e3;
volatile double BARYONIC_SCALE_LENGTH = 0; // fixed at main

volatile double SOFTENING_SPEED = 0;
volatile double SOFTENING_RADIUS = 0;

volatile double DARK_MATTER_SCALE_RADIUS = 0; // fixed at main
volatile double DARK_MATTER_DENSITY_0 = 0; // fixed at main

volatile double GAS_DENSITY = 0;
volatile double SIM_DT = 1e-6;

double darkMatterDensity0(double c)
{
    double factor = log(1 + c) - c / (1 + c);
    return HALO_MASS / (4 * PI * pow(DARK_MATTER_SCALE_RADIUS, 3) * factor);
}

double getNorm(double *vector)
{
    return pow(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2], 0.5);
}

double darkMatterDensity(double r)
{
    // r += SOFTENING_RADIUS;
    double factor = r / DARK_MATTER_SCALE_RADIUS;
    return DARK_MATTER_DENSITY_0 / (factor * pow(1 + factor, 2));
}

double darkMatterMass(double r)
{
    double factor = log(1 + r / DARK_MATTER_SCALE_RADIUS) - r / (DARK_MATTER_SCALE_RADIUS + r);
    return 4 * PI * DARK_MATTER_DENSITY_0 * factor * pow(DARK_MATTER_SCALE_RADIUS, 3);
}

double baryonicDensityHernquist(double r)
{
    // r += SOFTENING_RADIUS;
    return BARYONIC_TOTAL_MASS * BARYONIC_SCALE_LENGTH / (2 * PI * r * pow(r + BARYONIC_SCALE_LENGTH, 3));
}

double baryonicMassHernquist(double r)
{
    return BARYONIC_TOTAL_MASS * pow(r / (r + BARYONIC_SCALE_LENGTH), 2);
}

void sphericalToCartesian(double *r, double *theta, double *phi)
{
    double x, y, z;
    x = *r * sin(*theta) * cos(*phi);
    y = *r * sin(*theta) * sin(*phi);
    z = *r * cos(*theta);

    *r = x;
    *theta = y;
    *phi = z;
}

void printStatus(struct reb_simulation *sim, const char *filename)
{
    double t, x, y, z, vx, vy, vz;
    FILE *file;
    struct reb_particle particle = sim -> particles[0];
    t = sim -> t;
    x = particle.x;
    y = particle.y;
    z = particle.z;
    vx = particle.vx;
    vy = particle.vy;
    vz = particle.vz;
    file = fopen(filename, "a");
    fprintf(file, "%e %e %e %e %e %e %e %e\n", t, x, y, z, vx, vy, vz, SMBH_MASS);
    fclose(file);
}

struct reb_simulation* setupSimulation(double mass, double *position, double *speed,
                            void (*additional_force)(struct reb_simulation*),
                            int velocity_dependent)
{
    struct reb_simulation* sim = reb_create_simulation();
    sim->G = G0;
    sim->additional_forces = additional_force;
    sim->force_is_velocity_dependent = velocity_dependent;

    struct reb_particle particle = {0};
    particle.m = mass;
    particle.x = position[0];
    particle.y = position[1];
    particle.z = position[2];
    particle.vx = speed[0];
    particle.vy = speed[1];
    particle.vz = speed[2];
    reb_add(sim, particle);
    SMBH_MASS = mass;

    return sim;
}

void runSimulation(struct reb_simulation *sim, int n_points, int save_points, const char *filename)
{
    if (save_points > n_points) save_points = n_points;

    int i, every = n_points / save_points;

    sim->dt = SIM_DT;

    for (i = 0; i < n_points; i++)
    {
        if(i % every == 0)
        {
            printStatus(sim, filename);
        }
        // SIM_DT = sim->dt;
        reb_integrate(sim, sim->t + sim->dt);
    }
}

double getLocalSoundSpeed(double z)
{
    double factor1, factor2;
    factor1 = HALO_MASS / pow(1e2, 1./3);
    factor2 = (MATTER_DENSITY_PARAMETER * pow(h, 2) / 0.14);
    return 1.8 * pow(1 + z, 0.5) * factor1 * factor2;
}

double gasDensity(double r)
{
    return GAS_DENSITY * pow(r + SOFTENING_RADIUS, -2);
}

double *dynamicalFrictionDM(double *position, double *speed)
{
    double r = getNorm(position);
    double v = getNorm(speed);

    double mass = HALO_MASS;
    // double mass = darkMatterMass(r);
    double sigma = pow(0.5 * G0 * mass / R_VIR, 0.5);
    double x = v / (pow(2, 0.5) * sigma);

    double rho = darkMatterDensity(r) + baryonicDensityHernquist(r);
    double factor = -4 * PI * pow(G0, 2) * SMBH_MASS * rho * LN_LAMBDA
                * (erf(x) - 2 / pow(PI, 0.5) * x * exp(-pow(x, 2)));

    factor *= 1 / (pow(v, 3) + SOFTENING_RADIUS);
    // factor *= 0.5;
    double *ac = malloc(3 * sizeof(double));

    ac[0] = factor * speed[0];
    ac[1] = factor * speed[1];
    ac[2] = factor * speed[2];
    return ac;
}

double *dynamicalFrictionGas(double *position, double *speed)
{
    double f, all;
    double cs = getLocalSoundSpeed(20);
    double v = getNorm(speed);
    double mach = v / cs;
    if (mach <= 1.7)
    {
        double factor = erf(mach / pow(2, 0.5)) - pow(2 / PI, 0.5) * mach * exp(-0.5 * pow(mach, 2));
        if (mach <= 0.8) f = 0.5 * LN_LAMBDA * factor;
        else f = 1.5 * LN_LAMBDA * factor;
    }
    else f = 0.5 * log(1 - pow(mach, -2)) + LN_LAMBDA;

    double rho = gasDensity(getNorm(position));

    all = -4 * PI * pow(G0, 2) * SMBH_MASS * rho * f / (pow(v, 3) + SOFTENING_SPEED);
    double *ac = malloc(3 * sizeof(double));
    ac[0] = all * speed[0];
    ac[1] = all * speed[1];
    ac[2] = all * speed[2];
    return  ac;
}

double SMBHAccretion(double *position, double *speed)
{
    double r = getNorm(position);
    double v = pow(pow(getLocalSoundSpeed(20), 2) + pow(getNorm(speed), 2), 1.5);

    double bondi = 4 * PI * pow(G0 * SMBH_MASS, 2) * baryonicDensityHernquist(r) / v;
    double eddington = (1 - 0.1) * SMBH_MASS / (0.1 * 0.44);
    if (bondi < eddington) return bondi;
    return eddington;
}

double gravitationalForce(double r)
{
    double m = darkMatterMass(r);
    m += baryonicMassHernquist(r);
    return -G0 * m / pow(r, 2);
}

void baseCase(struct reb_simulation* sim)
{
    struct reb_particle *particle = &(sim->particles[0]);
    double pos[3] = {particle->x, particle->y, particle->z};
    double speed[3] = {particle->vx, particle->vy, particle->vz};
    double r = getNorm(pos);
    double dir_[3] = {pos[0] / r, pos[1] / r, pos[2] / r};

    double grav = gravitationalForce(r);
    double *df_g = dynamicalFrictionGas(pos, speed);
    double *df_dm = dynamicalFrictionDM(pos, speed);
    double df[3] = {df_g[0] + df_dm[0], df_g[1] + df_dm[1], df_g[2] + df_dm[2]};

    free(df_g);
    free(df_dm);

    double m_change = SMBHAccretion(pos, speed);
    SMBH_MASS += m_change * SIM_DT;

    double v = getNorm(speed);
    double accretion = v * m_change / SMBH_MASS;

    grav += accretion;

    particle->ax = grav * dir_[0] + df[0];
    particle->ay = grav * dir_[1] + df[1];
    particle->az = grav * dir_[2] + df[2];
}

void setR_vir(double r)
{
  R_VIR = r;
  // SOFTENING_SPEED = r * 1e-5;
  // SOFTENING_RADIUS = r * 1e-5;
  DARK_MATTER_SCALE_RADIUS = R_VIR / CONCENTRATION_PARAMETER;
  DARK_MATTER_DENSITY_0 = darkMatterDensity0(CONCENTRATION_PARAMETER);
  BARYONIC_SCALE_LENGTH = 0.01 * R_VIR / (1 + pow(2, 0.5));
}

void printConstants(void)
{
  printf("R_VIR: %f\nDARK_MATTER_SCALE_RADIUS: %f\nDARK_MATTER_DENSITY_0: %f\nBARYONIC_SCALE_LENGTH: %f\n\n",
            R_VIR, DARK_MATTER_SCALE_RADIUS,
            DARK_MATTER_DENSITY_0, BARYONIC_SCALE_LENGTH);
}

void run(double *positions, double *speeds, double smbh_mass, double dt, int n_points, int n_saves, const char *filename)
{
    SIM_DT = dt;
    struct reb_simulation* sim = setupSimulation(smbh_mass, positions, speeds, baseCase, 1);
    sim->integrator = REB_INTEGRATOR_LEAPFROG;
    runSimulation(sim, n_points, n_saves, filename);
}
