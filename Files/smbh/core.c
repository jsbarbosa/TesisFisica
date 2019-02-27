#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

volatile double Z = 20;

volatile double FB = 0.156;

volatile double GAS_CORE = 1e-3;
volatile double GAS_POWER = -2.2;
volatile double STELLAR_FRACTION = 0.0; // complement of gas percent

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

volatile double LAST_MAXIMA = -1;
volatile double LAST_SPEEDS[2] = {-1, -1};

volatile int STOP_SIMULATION = 0;

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
  STELLAR_FRACTION = ratio;
  setStellarTotalMass();
  setGasDensity();
}

void setStellarTotalMass(void)
{
  STELLAR_TOTAL_MASS = STELLAR_FRACTION * FB * HALO_MASS;
}

void setGasDensity(void)
{
  double m = 3 + GAS_POWER;
  double f1 = (pow(R_VIR, m) - pow(GAS_CORE, m)) / (m * pow(GAS_CORE, GAS_POWER));
  f1 = 4 * M_PI * (f1 + pow(GAS_CORE, 3) / 3);
  GAS_DENSITY = (1 - STELLAR_FRACTION) * FB * HALO_MASS / f1;
  // GAS_DENSITY =  / (4 * M_PI * pow(R_VIR, GAS_POWER + 3));
}

void setGasPower(double n)
{
  GAS_POWER = n;
  setGasDensity();
}

double darkMatterDensity0(double c)
{
  double factor = log(1 + c) - c / (1 + c);
  return DARK_MATTER_TOTAL_MASS / (4 * M_PI * pow(DARK_MATTER_SCALE_RADIUS, 3) * factor);
}

double darkMatterVelocityDispersion(void)
{
  double m = DARK_MATTER_TOTAL_MASS;
  return sqrt(0.5 * G0 * m / R_VIR);
}

double dampingFactor(double r, double v)
{
  double x = v / (sqrt(2) * darkMatterVelocityDispersion());
  x = erf(x) - (2 / sqrt(M_PI)) * x * exp(-x * x);
  return x;
}

double getNorm(double *vector)
{
  return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

double darkMatterDensity(double r)
{
  double factor = r / DARK_MATTER_SCALE_RADIUS;
  factor = DARK_MATTER_DENSITY_0 / (factor * pow(1 + factor, 2));
  if(factor > MAX_DENSITY_DM) return MAX_DENSITY_DM;
  return factor;
}

double darkMatterMass(double r)
{
  double factor = log(1 + r / DARK_MATTER_SCALE_RADIUS) - r / (DARK_MATTER_SCALE_RADIUS + r);
  return 4 * M_PI * DARK_MATTER_DENSITY_0 * factor * pow(DARK_MATTER_SCALE_RADIUS, 3);
}

double stellarDensityHernquist(double r)
{
  double factor = STELLAR_TOTAL_MASS * STELLAR_SCALE_LENGTH / (2 * M_PI * r * pow(r + STELLAR_SCALE_LENGTH, 3));
  if(factor > MAX_DENSITY_STARS) return MAX_DENSITY_STARS;
  return factor;
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

double getLocalSoundSpeed(double z)
{
  return SOUND_SPEED_FACTOR * sqrt(HALO_MASS / R_VIR);
    // double factor1, factor2, h;
    // h = getHubbleParameter(z);
    // factor1 = pow(HALO_MASS / 1e2, 1./3);
    // factor2 = pow(MATTER_DENSITY_PARAMETER * pow(h, 2) / 0.14, 1/6.);
    // return 1.8 * pow(1 + z, 0.5) * factor1 * factor2;
}

double gasDensity(double r)
{
  if (r < GAS_CORE) return GAS_DENSITY;
  return GAS_DENSITY * pow(GAS_CORE / r, -GAS_POWER);
  // double factor = GAS_DENSITY * pow(r, GAS_POWER);
  // if(factor > MAX_DENSITY_GAS) return MAX_DENSITY_GAS;
  // return factor;
}

double gasMass(double r)
{
  if (r < GAS_CORE)
  {
    return 4 * M_PI * GAS_DENSITY * pow(r, 3) / 3;
  }
  double m = 3 + GAS_POWER;
  double f1 = (pow(r, m) - pow(GAS_CORE, m)) / (m * pow(GAS_CORE, GAS_POWER));
  return 4 * M_PI * GAS_DENSITY * (f1 + pow(GAS_CORE, 3) / 3);
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
  factor *= -4 * M_PI * (G0 * G0) * SMBH_MASS * rho * LN_LAMBDA / (v * v);
  return factor;
}

double machFunction(double m)
{
  double factor = erf(m / sqrt(2)) -sqrt(2.0 / M_PI) * m * exp(-0.5 * m * m);
  if (m <= 0.8) return 0.5 * LN_LAMBDA * factor;
  else if(m <= 1.73100478) return 1.5 * LN_LAMBDA * factor;
  return 0.5 * log(1 - pow(m, -2.0)) + LN_LAMBDA;
}

double dynamicalFrictionGas(double r, double v)
{
    double cs = getLocalSoundSpeed(Z);
    double mach = v / cs;
    double f  = machFunction(mach);
    double rho = gasDensity(r);
    return -4 * M_PI * pow(G0, 2) * SMBH_MASS * rho * f / (v * v);
}

double SMBHAccretion(double *position, double *speed)
{
    double r = getNorm(position);
    double v = getNorm(speed);
    r = sqrt(r * r + SOFTENING_RADIUS);
    v = sqrt(v * v + SOFTENING_SPEED);
    v = pow(pow(getLocalSoundSpeed(Z), 2) + pow(v, 2), 1.5);
    double rho = gasDensity(r); // + stellarDensityHernquist(r)
    double bondi = 4 * M_PI * pow(G0 * SMBH_MASS, 2) * rho / v;
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

void localMaxima(double r, double v, double sim_time)
{
  if (LAST_SPEEDS[0] + LAST_SPEEDS[1] > 0)
  {
    if((LAST_SPEEDS[1] >= LAST_SPEEDS[0]) & (LAST_SPEEDS[1] >= v)) // is a local maxima
    {
      if ((v > LAST_MAXIMA) & (sim_time > 1e-3))
      {
        STOP_SIMULATION = 1;
        return ;
      }
      //
      // else if(r < 1e-3 * R_VIR)
      // {
      //   STOP_SIMULATION = 1;
      //   return ;
      // }
      LAST_MAXIMA = v;
    }
    LAST_SPEEDS[0] = LAST_SPEEDS[1];
    LAST_SPEEDS[1] = v;
  }
  else if(LAST_SPEEDS[0] == -1)
  {
    LAST_MAXIMA = v;
    LAST_SPEEDS[0] = v;
  }
  else
  {
    LAST_SPEEDS[1] = v;
  }
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

double getR_vir(void)
{
  return R_VIR;
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
  return pow((G * HALO_MASS) / (100 * pow(H, 2)), 1/3.);
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

    localMaxima(r, v, sim->t);
}

void printConstants(void)
{
  printf("R_VIR = %f\n", R_VIR);

  printf("\n##### FRACTIONS #####\n");
  printf("FB = %f\n", FB);
  printf("STELLAR_FRACTION = %f\n", STELLAR_FRACTION);

  printf("\n##### SCALES #####\n");
  printf("DARK_MATTER_SCALE_RADIUS = %f\n", DARK_MATTER_SCALE_RADIUS);
  printf("STELLAR_SCALE_LENGTH = %f\n", STELLAR_SCALE_LENGTH);

  printf("\n##### MASSES #####\n");
  printf("DARK_MATTER_TOTAL_MASS = %f\n", DARK_MATTER_TOTAL_MASS);
  printf("STELLAR_TOTAL_MASS = %f\n", STELLAR_TOTAL_MASS);

  printf("\n##### DENSITIES #####\n");
  printf("DARK_MATTER_DENSITY_0 = %f\n", DARK_MATTER_DENSITY_0);
  printf("GAS_DENSITY = %f\n", GAS_DENSITY);
  printf("GAS_POWER = %f\n", GAS_POWER);
}

void printStatus(struct reb_simulation* sim, const char *filename, int header)
{
  FILE *file;

  if (header == 0)
  {
    file = fopen(filename, "w");

    fprintf(file, "R_VIR = %f\n", R_VIR);
    fprintf(file, "FB = %f\n", FB);
    fprintf(file, "STELLAR_FRACTION = %f\n", STELLAR_FRACTION);

    fprintf(file, "DARK_MATTER_SCALE_RADIUS = %f\n", DARK_MATTER_SCALE_RADIUS);
    fprintf(file, "STELLAR_SCALE_LENGTH = %f\n", STELLAR_SCALE_LENGTH);

    fprintf(file, "DARK_MATTER_TOTAL_MASS = %f\n", DARK_MATTER_TOTAL_MASS);
    fprintf(file, "STELLAR_TOTAL_MASS = %f\n", STELLAR_TOTAL_MASS);

    fprintf(file, "DARK_MATTER_DENSITY_0 = %f\n", DARK_MATTER_DENSITY_0);
    fprintf(file, "GAS_DENSITY = %f\n", GAS_DENSITY);
    fprintf(file, "GAS_POWER = %f\n", GAS_POWER);

    // fprintf(file, "INTEGRATOR = %f\n", GAS_POWER);

    fprintf(file, "time\tz\tx\ty\tz\tvx\tvy\tvz\tmass\tlyapunov\n");
    fclose(file);
  }

  file = fopen(filename, "a");
  double t, x, y, z, vx, vy, vz, ly;
  struct reb_particle  particle = sim -> particles[0];
  t = sim -> t;
  x = particle.x;
  y = particle.y;
  z = particle.z;
  vx = particle.vx;
  vy = particle.vy;
  vz = particle.vz;
  ly = reb_tools_calculate_lyapunov(sim);
  fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", t, x, y, z, vx, vy, vz, SMBH_MASS, ly);
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

    reb_tools_megno_init(sim);
    return sim;
}

void runSimulation(struct reb_simulation* sim, int save_every, const char *filename)
{
  int i = 0;

  sim->dt = SIM_DT;
  STOP_SIMULATION = 0;
  LAST_MAXIMA = -1;
  LAST_SPEEDS[0] = -1;
  LAST_SPEEDS[1] = -1;

  while((!STOP_SIMULATION) & (sim->t + T0 < 13.78))
  {
    if(i % save_every == 0)
    {
      printStatus(sim, filename, i);
    }
    reb_integrate(sim, sim->t + sim->dt);
    i += 1;
  }
}

void run(double *positions, double *speeds, double smbh_mass, double dt, int integrator, int save_every, const char *filename)
{
  SIM_DT = dt;
  struct reb_simulation* sim = setupSimulation(smbh_mass, positions, speeds, baseCase);
  switch(integrator)
  {
    case INT_LEAPFROG:
      sim->integrator = REB_INTEGRATOR_LEAPFROG;
      break;
    case INT_IAS15:
      sim->integrator = REB_INTEGRATOR_IAS15;
      if(dt > 1e-10) SIM_DT = 1e-10;
      break;
    case INT_WHFAST:
      sim->integrator = REB_INTEGRATOR_WHFAST;
      break;
    case INT_JANUS:
      sim->integrator = REB_INTEGRATOR_JANUS;
      break;
    case INT_MERCURIUS:
      sim->integrator = REB_INTEGRATOR_MERCURIUS;
      break;
    default :
      sim->integrator = REB_INTEGRATOR_WHFAST;
 }

  runSimulation(sim, save_every, filename);
}
