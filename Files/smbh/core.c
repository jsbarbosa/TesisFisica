#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

volatile double Z = 20;

volatile double FB = 0.156;

volatile double GAS_CORE = 1e-3;
volatile double GAS_POWER = 2.2;
volatile double STELLAR_FRACTION = 0.01; // complement of gas percent

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
volatile double LAST_POSITIONS[2] = {-1, -1};

volatile unsigned int IS_LEAVING = 1;

volatile int STOP_SIMULATION = 0;

double const Z_TIME_COEFFS[Z_TIME_DEGREE + 1] = {-2.22289277, 5.13671586, -4.92900515, 3.71708597};
double const Z_HUBBLE_COEFFS[Z_HUBBLE_DEGREE + 1] = {0.0039385, 0.11201693, -0.11131262};

volatile double TRIAXIAL_A_1 = 1;
volatile double TRIAXIAL_A_2 = 0.5;
volatile double TRIAXIAL_A_3 = 0.25;

volatile double LYAPUNOV_DISTANCE_X = 1e-5;
volatile double LYAPUNOV_DISTANCE_P = 0;

volatile double RETURN_PROPERTIES[3] = {0, 0, 0};

const double ROOTS[GAUSS_DEGREE] = {-0.9988664 , -0.99403197, -0.98535408, -0.97286439, -0.95661096,
        -0.93665662, -0.91307856, -0.88596798, -0.85542977, -0.82158207,
        -0.78455583, -0.7444943 , -0.70155247, -0.65589647, -0.60770293,
        -0.5571583 , -0.50445814, -0.44980633, -0.39341431, -0.33550025,
        -0.27628819, -0.21600724, -0.15489059, -0.0931747 , -0.03109834,
         0.03109834,  0.0931747 ,  0.15489059,  0.21600724,  0.27628819,
         0.33550025,  0.39341431,  0.44980633,  0.50445814,  0.5571583 ,
         0.60770293,  0.65589647,  0.70155247,  0.7444943 ,  0.78455583,
         0.82158207,  0.85542977,  0.88596798,  0.91307856,  0.93665662,
         0.95661096,  0.97286439,  0.98535408,  0.99403197,  0.9988664 };
const double WEIGHTS[GAUSS_DEGREE] = {0.00290862, 0.0067598 , 0.01059055, 0.01438082, 0.01811556,
        0.02178024, 0.02536067, 0.02884299, 0.03221373, 0.03545984,
        0.03856876, 0.04152846, 0.0443275 , 0.04695505, 0.04940094,
        0.0516557 , 0.05371062, 0.05555774, 0.05718993, 0.05860085,
        0.05978506, 0.06073797, 0.0614559 , 0.06193607, 0.06217662,
        0.06217662, 0.06193607, 0.0614559 , 0.06073797, 0.05978506,
        0.05860085, 0.05718993, 0.05555774, 0.05371062, 0.0516557 ,
        0.04940094, 0.04695505, 0.0443275 , 0.04152846, 0.03856876,
        0.03545984, 0.03221373, 0.02884299, 0.02536067, 0.02178024,
        0.01811556, 0.01438082, 0.01059055, 0.0067598 , 0.00290862};

struct Results
{
  char **variables;
  double *values;
  double *t;
  double *x;
  double *y;
  double *z;
  double *vx;
  double *vy;
  double *vz;
  double *mass;

  int n_values, n_points;
};

double getG(void)
{
  return G0;
}

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
  // double m = 3 + GAS_POWER;
  // double f1 = (pow(R_VIR, m) - pow(GAS_CORE, m)) / (m * pow(GAS_CORE, GAS_POWER));
  // f1 = 4 * M_PI * (f1 + pow(GAS_CORE, 3) / 3);
  // GAS_DENSITY = (1 - STELLAR_FRACTION) * FB * HALO_MASS / f1;
  GAS_DENSITY = 1;
  GAS_DENSITY = (1 - STELLAR_FRACTION) * FB * HALO_MASS / gasMass(R_VIR);
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

double getM(double x, double y, double z)
{
  return sqrt(x * x + pow(TRIAXIAL_A_1 * y / TRIAXIAL_A_2, 2) + pow(TRIAXIAL_A_1 * z / TRIAXIAL_A_3, 2));
}

double getMTau(double x, double y, double z, double tau)
{
  x = x * x / (TRIAXIAL_A_1 * TRIAXIAL_A_1 + tau);
  y = y * y / (TRIAXIAL_A_2 * TRIAXIAL_A_2 + tau);
  z = z * z / (TRIAXIAL_A_3 * TRIAXIAL_A_3 + tau);
  return TRIAXIAL_A_1 * sqrt(x + y + z);
}

double darkMatterDensity(double r)
{
  double factor = r / DARK_MATTER_SCALE_RADIUS;
  factor = DARK_MATTER_DENSITY_0 / (factor * pow(1 + factor, 2));
  if(factor > MAX_DENSITY_DM) return MAX_DENSITY_DM;
  return factor;
}

double darkMatterDensityTriaxial(double x, double y, double z)
{
  double m = getM(x, y, z);
  return darkMatterDensity(m);
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
}

double gasDensity(double r)
{
  // if (r < GAS_CORE) return GAS_DENSITY;
  // return GAS_DENSITY * pow(GAS_CORE / r, -GAS_POWER);
  return GAS_DENSITY * pow((r / GAS_CORE + 1), - GAS_POWER);
}

double gasMass(double r)
{
  double u = r / GAS_CORE;
  double factor = -pow(GAS_POWER, 2) * (u + 1) * pow(u, 2) + 2 * (pow(u + 1, GAS_POWER) - pow(u, 3) - 1)
                + GAS_POWER * u * (3 * pow(u, 2) + u - 2);
  factor *= 4 * M_PI * GAS_DENSITY * pow(GAS_CORE, 3) * pow(u + 1, -GAS_POWER);
  return factor / ((GAS_POWER - 3) * (GAS_POWER - 2) * (GAS_POWER - 1));
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

double SMBHAccretion(double r, double v)
{
    v = pow(pow(getLocalSoundSpeed(Z), 2) + pow(v, 2), 1.5);
    double rho = gasDensity(r);
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

double *phiVector(double x, double y, double z, double tau)
{
  int i;
  double pos[3] = {x, y, z};
  double *phi = malloc(3 * sizeof(double));
  double a[3] = {TRIAXIAL_A_1, TRIAXIAL_A_2, TRIAXIAL_A_3};

  double factor = 1;
  for(i = 0; i < 3; i++) factor *= sqrt(tau + a[i] * a[i]);
  for(i = 0; i < 3; i++)
  {
    phi[i] = pos[i] / (factor * (tau + a[i] * a[i]));
  }
  return phi;
}

double *triaxial_gravitationalDarkMatter(double x, double y, double z, double tau)
{
  int i;
  double *phi = phiVector(x, y, z, tau);
  double m = getMTau(x, y, z, tau);

  double factor = m * pow(DARK_MATTER_SCALE_RADIUS + m, 2);
  for(i = 0; i < 3; i++) phi[i] *= 1.0 / factor;
  return phi;
}

double *triaxial_gravitationalStellar(double x, double y, double z, double tau)
{
  int i;
  double *phi = phiVector(x, y, z, tau);
  double m = getMTau(x, y, z, tau);

  double factor = m * pow(STELLAR_SCALE_LENGTH + m, 3);
  for(i = 0; i < 3; i++) phi[i] *= 1.0 / factor;
  return phi;
}

double *triaxial_gravitationalGas(double x, double y, double z, double tau)
{
  int i;
  double *phi = phiVector(x, y, z, tau);
  double m = getMTau(x, y, z, tau);
  for(i = 0; i < 3; i++)
  {
    phi[i] *= pow(GAS_CORE / (m + GAS_CORE), GAS_POWER);
  }
  return phi;
}

double *simpson(double *(*func)(double, double, double, double), double x, double y, double z, double gamma)
{
  int i, j, coeff;
  double value, h, tau, omega, diff;

  h = (1.0) / (INT_STEPS - 1.0);
  // h = 0.5 * M_PI / (INT_STEPS - 1.0);

  double *temp;
  double *grad = calloc(3, sizeof(double));
  for(i = 0; i < INT_STEPS; i++)
  {
    if ((i == 0) | (i == INT_STEPS - 1)) coeff = 1;
    else if(i % 3 == 0) coeff = 2;
    else coeff = 3;

    omega = h * i;
    tau = pow(omega / (1.0 - omega), 1.0 / gamma);
    diff = tau / (omega * gamma * (1 - omega));
    if(!isnormal(diff) & (diff != 0)) diff = 0;
    temp = func(x, y, z, tau);

    for(j = 0; j < 3; j++)
    {
      value = temp[j];
      if (isnormal(value) | (value == 0))
      {
        grad[j] += coeff * value * diff;
      }
    }
    free(temp);
  }

  for(i = 0; i < 3; i++) grad[i] *= (3 * h / 8) * 0.9988627714549503;
  return grad;
}

double *gaussLegendre(double *(*func)(double, double, double, double), double x, double y, double z, double gamma)
{
  int i, j;
  double tau, omega, diff, value;

  double *temp;
  double *grad = calloc(3, sizeof(double));
  for(i = 0; i < GAUSS_DEGREE; i++)
  {
    omega = 0.5 * ROOTS[i] + 0.5;
    tau = pow(omega / (1.0 - omega), 1.0 / gamma);
    diff = tau / (omega * gamma * (1 - omega));
    if(!isnormal(diff) & (diff != 0)) diff = 0;
    temp = func(x, y, z, tau);

    for(j = 0; j < 3; j++)
    {
      value = temp[j];
      if (isnormal(value) | (value == 0))
      {
        grad[j] += value * diff * WEIGHTS[i];
      }
    }
    free(temp);
  }

  for(i = 0; i < 3; i++) grad[i] *= 0.5;
  return grad;
}

double *triaxial_gravDM(double x, double y, double z, double gamma)
{
  int i;
  double *grad = gaussLegendre(triaxial_gravitationalDarkMatter, x, y, z, gamma);
  for(i = 0; i < 3; i++) grad[i] *= 2 * M_PI * G0 * pow(DARK_MATTER_SCALE_RADIUS, 3) *
          DARK_MATTER_DENSITY_0 * TRIAXIAL_A_1 * TRIAXIAL_A_2 * TRIAXIAL_A_3;
  return grad;
}

double *triaxial_gravS(double x, double y, double z, double gamma)
{
  int i;
  double *grad = gaussLegendre(triaxial_gravitationalStellar, x, y, z, gamma);
  for(i = 0; i < 3; i++) grad[i] *= G0 * STELLAR_TOTAL_MASS * TRIAXIAL_A_1 * TRIAXIAL_A_2 * TRIAXIAL_A_3 * STELLAR_SCALE_LENGTH;
  return grad;
}

double *triaxial_gravG(double x, double y, double z, double gamma)
{
  int i;
  double *grad = gaussLegendre(triaxial_gravitationalGas, x, y, z, gamma);
  for(i = 0; i < 3; i++) grad[i] *= 2 * M_PI * G0 * GAS_DENSITY * TRIAXIAL_A_1 * TRIAXIAL_A_2 * TRIAXIAL_A_3;
  return grad;
}

int localMaxima(double r, double sim_time)
{
  int value = 0;
  // if((r <= RETURN_FRACTION * R_VIR) & (LAST_MAXIMA < 2 * RETURN_FRACTION * R_VIR) & (RETURN_PROPERTIES[0] < 0))
  double diff = fabs(r / R_VIR - RETURN_FRACTION);
  if(diff < 1e-4)
  {
    // printf("%f, %f, %f\n", r / R_VIR, RETURN_FRACTION, diff);
    value = 1;
  }
  if (LAST_POSITIONS[0] + LAST_POSITIONS[1] > 0)
  {
    if((LAST_POSITIONS[1] >= LAST_POSITIONS[0]) & (LAST_POSITIONS[1] >= r)) // is a local maxima
    {
      if(LAST_MAXIMA > 0)
      {
        if ((r > 1.25 * LAST_MAXIMA) & (sim_time > 1e-3))
        {
          STOP_SIMULATION = 1;
          return value;
        }
      }
      LAST_MAXIMA = r;
    }
    if(LAST_MAXIMA > 0)
    {
      if(LAST_MAXIMA < 1e-3 * R_VIR)
      {
        STOP_SIMULATION = 1;
        return value;
      }
    }

    LAST_POSITIONS[0] = LAST_POSITIONS[1];
    LAST_POSITIONS[1] = r;
  }
  else if(LAST_POSITIONS[0] == -1)
  {
    LAST_POSITIONS[0] = r;
  }
  else
  {
    LAST_POSITIONS[1] = r;
  }

  return value;
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

void setTriaxialCoeffs(double x, double y, double z)
{
  TRIAXIAL_A_1 = x;
  TRIAXIAL_A_2 = y;
  TRIAXIAL_A_3 = z;
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

void sphericalCase(struct reb_simulation* sim)
{
    Z = getRedshift(sim->t + T0);

    struct reb_particle* particle = &(sim->particles[0]);
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

    double m_change = SMBHAccretion(r, v);
    SMBH_MASS += m_change * SIM_DT;
    particle->m = SMBH_MASS;

    double accretion = - v * m_change / SMBH_MASS;

    ax += (accretion + df) * dir_v[0];
    ay += (accretion + df) * dir_v[1];
    az += (accretion + df) * dir_v[2];

    particle->ax = ax;
    particle->ay = ay;
    particle->az = az;

    // sim->ri_whfast.recalculate_coordinates_this_timestep = 1;

    if(localMaxima(r, sim->t))
    {
      RETURN_PROPERTIES[0] = sim->t;
      RETURN_PROPERTIES[1] = particle->m;
    }
}

void triaxialCase(struct reb_simulation* sim)
{
    Z = getRedshift(sim->t + T0);

    struct reb_particle* particle = &(sim->particles[0]);
    double pos[3] = {particle->x, particle->y, particle->z};
    double speed[3] = {particle->vx, particle->vy, particle->vz};

    double ax, ay, az;
    double r = getNorm(pos);
    double m = getSoftenedLength(getM(pos[0], pos[1], pos[2]));
    double v = getSoftenedSpeed(getNorm(speed));
    double dir_v[3] = {speed[0] / v, speed[1] / v, speed[2] / v};

    double *p_dm = triaxial_gravDM(pos[0], pos[1], pos[2], 0.2);
    double *p_stars = triaxial_gravS(pos[0], pos[1], pos[2], 0.2);
    double *p_gas = triaxial_gravG(pos[0], pos[1], pos[2], 0.2);

    int i;
    for(i = 0; i < 3; i++) p_dm[i] += p_stars[i] + p_gas[i];

    ax = - p_dm[0];
    ay = - p_dm[1];
    az = - p_dm[2];

    free(p_dm);
    free(p_stars);
    free(p_gas);

    double df_dm = dynamicalFrictionDM(m, v);
    double df_g = dynamicalFrictionGas(m, v);
    double df = df_dm + df_g;

    double m_change = SMBHAccretion(m, v);
    SMBH_MASS += m_change * SIM_DT;

    double accretion = - v * m_change / SMBH_MASS;

    ax += (accretion + df) * dir_v[0];
    ay += (accretion + df) * dir_v[1];
    az += (accretion + df) * dir_v[2];

    particle->ax = ax;
    particle->ay = ay;
    particle->az = az;

    // sim->ri_whfast.recalculate_coordinates_this_timestep = 1;

    if(localMaxima(r, sim->t))
    {
      RETURN_PROPERTIES[0] = sim->t;
      RETURN_PROPERTIES[1] = particle->m;
    }
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

    fprintf(file, "time\tz\tx\ty\tz\tvx\tvy\tvz\tmass\n");
    fclose(file);
  }

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
  fprintf(file, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", t, x, y, z, vx, vy, vz, SMBH_MASS);
  fclose(file);
}

struct reb_simulation* setupSimulation(double mass, double *position, double *speed, int integrator,
                            void (*additional_force)(struct reb_simulation * ))
{
  int coeff = 1;
  struct reb_simulation* sim = reb_create_simulation();
  sim->G = G0;
  sim->gravity = REB_GRAVITY_NONE;
  sim->additional_forces = additional_force;
  sim->force_is_velocity_dependent = 1;

  struct reb_particle particle = {0};
  particle.m = mass;
  particle.x = coeff * position[0];
  particle.y = coeff * position[1];
  particle.z = coeff * position[2];
  particle.vx = coeff * speed[0];
  particle.vy = coeff * speed[1];
  particle.vz = coeff * speed[2];

  reb_add(sim, particle);

  if(coeff == 2)
  {
    struct reb_particle particle2 = {0};
    particle2.m = mass;
    particle2.x = 0;
    particle2.y = 0;
    particle2.z = 0;
    particle2.vx = 0;
    particle2.vy = 0;
    particle2.vz = 0;

    reb_add(sim, particle2);
    reb_move_to_com(sim);
  }

  SMBH_MASS = mass;

  // reb_tools_megno_init(sim);
  sim->exact_finish_time = 0;

  switch(integrator)
  {
    case INT_LEAPFROG:
      sim->integrator = REB_INTEGRATOR_LEAPFROG;
      break;
    case INT_IAS15:
      sim->integrator = REB_INTEGRATOR_IAS15;
      break;
    case INT_WHFAST:
      SIM_DT *= 1. / 10;
      sim->integrator = REB_INTEGRATOR_WHFAST;
      break;
    case INT_JANUS:
      sim->integrator = REB_INTEGRATOR_JANUS;
      break;
    case INT_MERCURIUS:
      sim->integrator = REB_INTEGRATOR_MERCURIUS;
      break;
    default :
      sim->integrator = REB_INTEGRATOR_LEAPFROG;
  }

  return sim;
}

void reWriteFile(const char *filename)
{
  int length = 250;
  char new_name[length];
  strcpy(new_name, filename);

  new_name[strlen(new_name)] = 't';

  FILE *file_in = fopen(filename, "r");
  FILE *file_out = fopen(new_name, "w");
  char line_buffer[length]; // prepares a list of length chars to store a single line of the document
  char *line;

  fprintf(file_out, "RETURN_TIME = %e\nRETURN_MASS = %e\n", RETURN_PROPERTIES[0],
            RETURN_PROPERTIES[1]);

  while(fgets(line_buffer, length, file_in) != NULL) // reads up to length characters of the dataFile and stores them at the line_buffer
  {
    line = strtok(line_buffer, "\n");
    fprintf(file_out, "%s\n", line);
  }

  fclose(file_in);
  fclose(file_out);

  rename(new_name, filename);
  }

void runSimulation(struct reb_simulation* sim, int save_every, const char *filename)
{
  int i = 0;

  sim->dt = SIM_DT;
  STOP_SIMULATION = 0;
  LAST_MAXIMA = -1;
  LAST_POSITIONS[0] = -1;
  LAST_POSITIONS[1] = -1;
  RETURN_PROPERTIES[0] = -1;
  RETURN_PROPERTIES[1] = -1;

  // while(!STOP_SIMULATION)
  while((!STOP_SIMULATION) & (sim->t + T0 < 13.799))
  {
    if(i % save_every == 0)
    {
      printStatus(sim, filename, i);
    }
    reb_integrate(sim, sim->t + sim->dt);
    i += 1;
  }
  reb_free_simulation(sim);

  reWriteFile(filename);
}

void run(double *positions, double *speeds, double smbh_mass, double dt, int triaxial, int integrator, int save_every, const char *filename)
{
  SIM_DT = dt;
  struct reb_simulation* sim;
  if(triaxial == 1) sim = setupSimulation(smbh_mass, positions, speeds, integrator, triaxialCase);
  else sim = setupSimulation(smbh_mass, positions, speeds, integrator, sphericalCase);

  runSimulation(sim, save_every, filename);
}

void printResults(struct Results results)
{
  int i;
  for(i = 0; i < results.n_values; i ++)
  {
    printf("%s = %f\n", results.variables[i], results.values[i]);
  }
  for(i = 0; i < results.n_points; i++)
  {
    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", results.t[i], results.x[i], results.y[i], results.z[i],
              results.vx[i], results.vy[i], results.vz[i], results.mass[i]);
  }
}

struct Results loadResults(const char *file_name)
{
  FILE *file;
  int length = 250;
  struct Results results;
  int i = 0, j = 0, k = 0;
  file = fopen(file_name, "r");

  char line_buffer[length]; // prepares a list of length chars to store a single line of the document
  char *split_buffer; // prepares a pointer of chars to store a single word or number in the line_buffer
  char *temp_buffer;

  double value;

  if (file == NULL)
  {
    // verifies the file exists
    printf("Error Reading File\n");
    exit(0);
  }

  while(fgets(line_buffer, length, file) != NULL) // reads up to length characters of the dataFile and stores them at the line_buffer
  {
    j = 0;
    split_buffer = strtok(line_buffer, "\t"); // first word of the line
    while (split_buffer != NULL)
    {
        temp_buffer = strtok(NULL, "\t");
        if((j == 0) & (temp_buffer == NULL))
        {
          split_buffer = strtok(split_buffer, " = ");
          temp_buffer = strtok(NULL, " = ");

          if(i == 0)
          {
            results.variables = calloc(1, sizeof(char *));
            results.values = calloc(1, sizeof(double));
          }
          else
          {
            results.variables = realloc(results.variables, sizeof(char *) * (i + 1));
            results.values = realloc(results.values, sizeof(double) * (i + 1));
          }

          results.variables[i] = calloc(50, sizeof(char));
          results.values[i] = atof(temp_buffer);
          strcpy(results.variables[i],  split_buffer);
          temp_buffer = NULL;

          results.n_values = i;
        }
        else
        {
          if (k == 0)
          {
            k += 1;
            break;
          }
          else
          {
            value = atof(split_buffer);
            switch (j)
            {
              case 0:
              {
                if(k == 2)
                {
                  results.t = calloc(1, sizeof(double));
                  results.x = calloc(1, sizeof(double));
                  results.y = calloc(1, sizeof(double));
                  results.z = calloc(1, sizeof(double));
                  results.vx = calloc(1, sizeof(double));
                  results.vy = calloc(1, sizeof(double));
                  results.vz = calloc(1, sizeof(double));
                  results.mass = calloc(1, sizeof(double));
                }
                else
                {
                  results.t = realloc(results.t, (k - 1) * sizeof(double));
                  results.x = realloc(results.x, (k - 1) * sizeof(double));
                  results.y = realloc(results.y, (k - 1) * sizeof(double));
                  results.z = realloc(results.z, (k - 1) * sizeof(double));
                  results.vx = realloc(results.vx, (k - 1) * sizeof(double));
                  results.vy = realloc(results.vy, (k - 1) * sizeof(double));
                  results.vz = realloc(results.vz, (k - 1) * sizeof(double));
                  results.mass = realloc(results.mass, (k - 1) * sizeof(double));
                }
                results.t[k - 2] = value;
                break;
              }
              case 1:
              {
                results.x[k - 2] = value;
                break;
              }
              case 2:
              {
                results.y[k - 2] = value;
                break;
              }
              case 3:
              {
                results.z[k - 2] = value;
                break;
              }
              case 4:
              {
                results.vx[k - 2] = value;
                break;
              }
              case 5:
              {
                results.vy[k - 2] = value;
                break;
              }
              case 6:
              {
                results.vz[k - 2] = value;
                break;
              }
              case 7:
              {
                results.mass[k - 2] = value;
                break;
              }
              default:
              {

              }
            }
          }
        }

        split_buffer = temp_buffer;
        j += 1;
    }
    if(k > 0) k += 1;
    i += 1;
  }

  results.n_points = k - 2;

  for(i = 0; i < results.n_values; i++)
  {
    value = results.values[i];
    if(strcmp(results.variables[i], "R_VIR") == 0) R_VIR = value;
    else if(strcmp(results.variables[i], "FB") == 0) FB = value;
    else if(strcmp(results.variables[i], "STELLAR_FRACTION") == 0) STELLAR_FRACTION = value;
    else if(strcmp(results.variables[i], "DARK_MATTER_SCALE_RADIUS") == 0) DARK_MATTER_SCALE_RADIUS = value;
    else if(strcmp(results.variables[i], "STELLAR_SCALE_LENGTH") == 0) STELLAR_SCALE_LENGTH = value;
    else if(strcmp(results.variables[i], "DARK_MATTER_TOTAL_MASS") == 0) DARK_MATTER_TOTAL_MASS = value;
    else if(strcmp(results.variables[i], "STELLAR_TOTAL_MASS") == 0) STELLAR_TOTAL_MASS = value;
    else if(strcmp(results.variables[i], "DARK_MATTER_DENSITY_0") == 0) DARK_MATTER_DENSITY_0 = value;
    else if(strcmp(results.variables[i], "GAS_DENSITY") == 0) GAS_DENSITY = value;
  }

  return results;
}

void getDelta(double **vals, double *refs)
{
  int i;
  for(i = 0; i < 3; i++)
  {
    (*vals)[i] = (*vals)[i] - refs[i];
  }
}

double getS(double *q, double *p, double *q_0, double *p_0)
{
  int i;
  double num = 0, dem = 0;
  for(i = 0; i < 3; i++)
  {
    num += pow(q[i], 2) + pow(p[i], 2);
    dem += pow(p_0[i], 2) + pow(q_0[i], 2);
  }

  return sqrt(num / dem);
}

void printTriplet(double *array)
{
  int i;
  for(i = 0; i < 3; i++)
  {
    printf("%e\t", array[i]);
  }
  printf("\n");
}

double lyapunov(double *positions, double *speeds, double *d_q0, double *d_p0,
        double smbh_mass, double T, double dt, int l, int triaxial)
{
  int i, j;
  struct reb_simulation* sim;
  struct reb_particle* particle;

  double ln = 0, s;
  double *q = malloc(3 * sizeof(double));
  double *p = malloc(3 * sizeof(double));
  double *ref_p = malloc(3 * sizeof(double));
  double *ref_q = malloc(3 * sizeof(double));

  for(i = 0; i < l + 1; i++)
  {
    if(i == 0)
    {
      if(triaxial == 1) sim = setupSimulation(smbh_mass, positions, speeds, INT_LEAPFROG, triaxialCase);
      else sim = setupSimulation(smbh_mass, positions, speeds, INT_LEAPFROG, sphericalCase);
    }
    else
    {
      for(j = 0; j < 3; j++)
      {
        q[j] = positions[j] + d_q0[j];
        // p[j] = speeds[j] + d_p0[j] / smbh_mass;
      }

      if(triaxial == 1) sim = setupSimulation(smbh_mass, q, speeds, INT_LEAPFROG, triaxialCase);
      else sim = setupSimulation(smbh_mass, q, speeds, INT_LEAPFROG, sphericalCase);
    }

    sim->dt = dt;
    reb_integrate(sim, T);

    particle = &(sim->particles[0]);
    q[0] = particle->x;
    q[1] = particle->y;
    q[2] = particle->z;
    p[0] = particle->vx * particle->m;
    p[1] = particle->vy * particle->m;
    p[2] = particle->vz * particle->m;
    smbh_mass = particle->m;

    reb_free_simulation(sim);

    if(i == 0)
    {
      for(j = 0; j < 3; j++)
      {
        ref_q[j] = q[j];
        ref_p[j] = p[j];
      }
      // printTriplet(q);
      // printTriplet(p);
    }

    else
    {
      // printTriplet(q);
      // printTriplet(p);
      getDelta(&q, ref_q);
      getDelta(&p, ref_p);

      s = getS(q, p, d_q0, d_p0);

      // printTriplet(q);
      // printTriplet(p);
      // printTriplet(d_q0);
      // printTriplet(d_p0);
      // printf("%d %f %f\n", i, s, log(s));
      ln += log(s);

      for(j = 0; j < 3; j++)
      {
        d_q0[j] = q[j] / s;
        // d_p0[j] = p[j] / s;
      }
    }
  }

  free(q);
  free(p);
  free(ref_p);
  free(ref_q);

  return ln / (l * T);
}

void testLoad(const char *file_name)
{
  struct Results results = loadResults(file_name);
  printResults(results);
}

double getReturnFraction(void)
{
  return RETURN_FRACTION;
}
