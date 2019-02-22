#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

volatile double Z = 20;

volatile double FB = 0.156;
volatile double STELLAR_RATIO = 0.0; // complement of gas percent

/* FIXED AT MAIN */
volatile double R_VIR = 0;
volatile double SMBH_MASS = 0;

volatile double DARK_MATTER_TOTAL_MASS = 0;
volatile double STELLAR_TOTAL_MASS = 0;
volatile double STELLAR_SCALE_LENGTH = 0;

volatile double SOFTENING_SPEED = 0;//1e-8;
volatile double SOFTENING_RADIUS = 1e-7;

volatile double DARK_MATTER_SCALE_RADIUS = 0; // fixed at main
volatile double DARK_MATTER_DENSITY_0 = 0; // fixed at main

volatile double GAS_DENSITY = 0; //fixed at main
volatile double SIM_DT;

double const Z_TIME_COEFFS[Z_TIME_DEGREE + 1] = {-2.22289277, 5.13671586, -4.92900515, 3.71708597};
double const Z_HUBBLE_COEFFS[Z_HUBBLE_DEGREE + 1] = {0.0039385, 0.11201693, -0.11131262};

void setBaryonicFraction(double fb)
{
  FB = fb;
  DARK_MATTER_TOTAL_MASS = (1 - FB) * HALO_MASS;
  DARK_MATTER_DENSITY_0 = darkMatterDensity0(CONCENTRATION_PARAMETER);
  setStellarTotalMass();
  setGasDensity();
}

void setStellarRatio(double ratio)
{
  STELLAR_RATIO = ratio;
  setStellarTotalMass();
  setGasDensity();
}

void setStellarTotalMass(void)
{
  STELLAR_TOTAL_MASS = STELLAR_RATIO * FB * HALO_MASS;
}

void setGasDensity(void)
{
  GAS_DENSITY = (1 - STELLAR_RATIO) * FB * HALO_MASS / (4 * PI * R_VIR);
}

double darkMatterDensity0(double c)
{
  double factor = log(1 + c) - c / (1 + c);
  return DARK_MATTER_TOTAL_MASS / (4 * PI * pow(DARK_MATTER_SCALE_RADIUS, 3) * factor);
}

double darkMatterVelocityDispersion(double r)
{
  // double m = darkMatterMass(r);
  double m = DARK_MATTER_TOTAL_MASS;
  return sqrt(0.5 * G0 * m / R_VIR);
}

double dampingFactor(double r, double v)
{
  double x = v / (sqrt(2) * darkMatterVelocityDispersion(r));
  x = erf(x) - (2 / sqrt(PI)) * x * exp(-x * x);
  return x;
}

double getNorm(double *vector)
{
  return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

double darkMatterDensity(double r)
{
  double factor = r / DARK_MATTER_SCALE_RADIUS;
  return DARK_MATTER_DENSITY_0 / (factor * pow(1 + factor, 2));
}

double darkMatterMass(double r)
{
  double factor = log(1 + r / DARK_MATTER_SCALE_RADIUS) - r / (DARK_MATTER_SCALE_RADIUS + r);
  return 4 * PI * DARK_MATTER_DENSITY_0 * factor * pow(DARK_MATTER_SCALE_RADIUS, 3);
}

double stellarDensityHernquist(double r)
{
  return STELLAR_TOTAL_MASS * STELLAR_SCALE_LENGTH / (2 * PI * r * pow(r + STELLAR_SCALE_LENGTH, 3));
}

double stellarMassHernquist(double r)
{
    return STELLAR_TOTAL_MASS * pow(r / (r + STELLAR_SCALE_LENGTH), 2);
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

void printStatus(struct reb_simulation* sim, const char *filename, int header)
{
  FILE *file;

  if (header == 0)
  {
    file = fopen(filename, "w");
    fprintf(file, "time\tz\tx\ty\tz\tvx\tvy\tvz\tmass\n");
  }
  else
  {
    file = fopen(filename, "a");
    double t, x, y, z, vx, vy, vz;
    struct reb_particle  particle = sim -> particles[0];
    t = sim -> t;
    x = particle.x;
    y = particle.y;
    z = particle.z;
    vx = particle.vx;
    vy = particle.vy;
    vz = particle.vz;
    fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", t, Z, x, y, z, vx, vy, vz, SMBH_MASS);
  }
  fclose(file);
}

struct reb_simulation* setupSimulation(double mass, double *position, double *speed,
                            void (*additional_force)(struct reb_simulation * ))
{
    struct reb_simulation* sim = reb_create_simulation();
    sim->G = G0;
    sim->gravity = REB_GRAVITY_NONE;
    sim->additional_forces = additional_force;
    sim->force_is_velocity_dependent = 1;

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

void runSimulation(struct reb_simulation* sim, int n_points, int save_points, const char *filename)
{
    if (save_points > n_points) save_points = n_points;

    int i, every = n_points / save_points;

    sim->dt = SIM_DT;

    for (i = 0; i < n_points; i++)
    {
        if(i % every == 0)
        {
            printStatus(sim, filename, i);
        }
        reb_integrate(sim, sim->t + sim->dt);
    }
}

double getLocalSoundSpeed(double z)
{
    double factor1, factor2, h;
    h = getHubbleParameter(z);
    factor1 = pow(HALO_MASS / 1e2, 1./3);
    factor2 = pow(MATTER_DENSITY_PARAMETER * pow(h, 2) / 0.14, 1/6.);
    return 1.8 * pow(1 + z, 0.5) * factor1 * factor2;
}

double gasDensity(double r)
{
    return GAS_DENSITY * pow(r, -2);
}

double gasMass(double r)
{
  return 4 * PI * GAS_DENSITY * r;
}

double getSoftenedLength(double r)
{
  return r = sqrt(r * r + SOFTENING_RADIUS);
}

double getSoftenedSpeed(double v)
{
  return v = sqrt(v * v + SOFTENING_SPEED);
}

double dynamicalFrictionDM(double r, double v)
{
  double factor = dampingFactor(r, v);
  double rho = darkMatterDensity(r) + stellarDensityHernquist(r);
  factor *= -4 * PI * (G0 * G0) * SMBH_MASS * rho * LN_LAMBDA / (v * v);
  return factor;
}

double dynamicalFrictionGas(double r, double v)
{
    double f;
    double cs = getLocalSoundSpeed(Z);
    double mach = v / cs;
    if (mach <= 1.7)
    {
      double factor = erf(mach / sqrt(2)) - sqrt(2 / PI) * mach * exp(-0.5 * mach * mach);
      if (mach <= 0.8) f = 0.5 * LN_LAMBDA * factor;
      else f = 1.5 * LN_LAMBDA * factor;
    }
    else f = 0.5 * log(1 - pow(mach, -2)) + LN_LAMBDA;

    double rho = gasDensity(r);
    return -4 * PI * pow(G0, 2) * SMBH_MASS * rho * f / (v * v * v);
}

double SMBHAccretion(double *position, double *speed)
{
    double r = getNorm(position);
    r = sqrt(r * r + SOFTENING_RADIUS);
    double v = getNorm(speed);
    v = sqrt(v * v + SOFTENING_SPEED);
    v = pow(pow(getLocalSoundSpeed(Z), 2) + pow(v, 2), 1.5);
    double rho = stellarDensityHernquist(r) + gasDensity(r);
    double bondi = 4 * PI * pow(G0 * SMBH_MASS, 2) * rho / v;
    double eddington = (1 - 0.1) * SMBH_MASS / (0.1 * 0.44);
    if (bondi < eddington) return bondi;
    return eddington;
}

double gravitationalForce(double r)
{
  double m = darkMatterMass(r);
  m += stellarMassHernquist(r);
  m += gasMass(r);
  return -G0 * m / pow(r, 2);
}

void baseCase(struct reb_simulation* sim)
{
    Z = getRedshift(sim->t + T0);

    struct reb_particle*  particle = &(sim->particles[0]);
    double pos[3] = {particle->x, particle->y, particle->z};
    double speed[3] = {particle->vx, particle->vy, particle->vz};
    double r = getSoftenedLength(getNorm(pos));
    double v = getSoftenedSpeed(getNorm(speed));

    double ax, ay, az;
    double dir_[3] = {pos[0] / r, pos[1] / r, pos[2] / r};
    double dir_v[3] = {speed[0] / v, speed[1] / v, speed[2] / v};

    double grav = gravitationalForce(r);

    ax = grav * dir_[0];
    ay = grav * dir_[1];
    az = grav * dir_[2];

    double df_dm = dynamicalFrictionDM(r, v);
    double df_g = dynamicalFrictionGas(r, v);
    double df = df_dm + df_g;

    double m_change = SMBHAccretion(pos, speed);
    SMBH_MASS += m_change * SIM_DT;

    double accretion = - v * m_change / SMBH_MASS;

    ax += (accretion + df) * dir_v[0];
    ay += (accretion + df) * dir_v[1];
    az += (accretion + df) * dir_v[2];

    particle->ax = ax;
    particle->ay = ay;
    particle->az = az;
}

void setR_vir(double r)
{
  R_VIR = r;

  setBaryonicFraction(FB);
  DARK_MATTER_SCALE_RADIUS = R_VIR / CONCENTRATION_PARAMETER;
  DARK_MATTER_DENSITY_0 = darkMatterDensity0(CONCENTRATION_PARAMETER);
  STELLAR_SCALE_LENGTH = 0.01 * R_VIR / (1 + pow(2, 0.5));

  setGasDensity();
}

void printConstants(void)
{
  printf("R_VIR: %f\n", R_VIR);

  printf("\n===== FRACTIONS =====\n");
  printf("FB: %f\n", FB);
  printf("STELLAR_RATIO: %f\n", STELLAR_RATIO);

  printf("\n===== SCALES =====\n");
  printf("DARK_MATTER_SCALE_RADIUS: %f\n", DARK_MATTER_SCALE_RADIUS);
  printf("STELLAR_SCALE_LENGTH: %f\n", STELLAR_SCALE_LENGTH);

  printf("\n===== MASSES =====\n");
  printf("DARK_MATTER_TOTAL_MASS: %f\n", DARK_MATTER_TOTAL_MASS);
  printf("STELLAR_TOTAL_MASS: %f\n", STELLAR_TOTAL_MASS);


  printf("\n===== DENSITIES =====\n");
  printf("DARK_MATTER_DENSITY_0: %f\n", DARK_MATTER_DENSITY_0);
  printf("GAS_DENSITY: %f\n", GAS_DENSITY);
}

void run(double *positions, double *speeds, double smbh_mass, double dt, int n_points, int n_saves, const char *filename)
{
    SIM_DT = dt;
    struct reb_simulation* sim = setupSimulation(smbh_mass, positions, speeds, baseCase);
    sim->integrator = REB_INTEGRATOR_LEAPFROG;
    runSimulation(sim, n_points, n_saves, filename);
}

double getRedshift(double t)
{
  int i;
  double z = 0;
  for(i = 0; i < Z_TIME_DEGREE + 1; i++)
  {
    z += Z_TIME_COEFFS[i] * pow(t, Z_TIME_DEGREE - i);
  }
  z = exp(z);
  return z;
}

double getHubbleParameter(double z)
{
  int i;
  double h = 0;
  for(i = 0; i < Z_HUBBLE_DEGREE + 1; i++)
  {
    h += Z_HUBBLE_COEFFS[i] * pow(z, Z_HUBBLE_DEGREE - i);
  }
  return h;
}

double calculateR_vir(double G, double H)
{
  return pow((G * 1e3) / (100 * pow(H, 2)), 1/3.);
}
