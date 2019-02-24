#define LN_LAMBDA 2.3

#define CONCENTRATION_PARAMETER 4
#define MATTER_DENSITY_PARAMETER 0.309 // LAMBDA_CDM model

#define SOUND_SPEED_FACTOR 0.6141441129704023
#define HALO_MASS 1e3

#define MAX_DENSITY_DM 1e8
#define MAX_DENSITY_STARS 1e8
#define MAX_DENSITY_GAS 1e8

#define G0 0.449849
// #define G0 4.49849e-06 // kpc gyr m0

#define Z_TIME_DEGREE 3
#define Z_HUBBLE_DEGREE 2
#define T0 0.18004611427373277

#define INT_LEAPFROG 0
#define INT_IAS15 1
#define INT_WHFAST 2
#define INT_SEI 3
#define INT_JANUS 4
#define INT_MERCURIUS 5

double getR_vir(void);
void setR_vir(double r);
void printConstants(void);
double gasDensity(double r);
double gasMass(double r);
double getNorm(double *vector);
double darkMatterMass(double r);
double machFunction(double mach);
double darkMatterDensity(double r);
double darkMatterDensity0(double c);
double getLocalSoundSpeed(double z);
double gravitationalForce(double r);
double baryonicMassHernquist(double r);
double baryonicDensityHernquist(double r);
void baseCase(struct reb_simulation* sim);
double SMBHAccretion(double *position, double *speed);
double dynamicalFrictionDM(double r, double v);
double dynamicalFrictionGas(double r, double v);

void sphericalToCartesian(double *r, double *theta, double *phi);
void printStatus(struct reb_simulation *sim, const char *filename, int header);

struct reb_simulation* setupSimulation(double mass, double *position, double *speed,
                            void (*additional_force)(struct reb_simulation*));
void runSimulation(struct reb_simulation *sim, int save_every, const char *filename);
void run(double *positions, double *speeds, double smbh_mass, double dt, int integrator, int save_every, const char *filename);

// void integrate(struct reb_simulation* sim);
// double *baseCase(double *position, double *speed);

double getRedshift(double t);
double getHubbleParameter(double z);
double calculateR_vir(double G, double H);
double darkMatterVelocityDispersion(void);
void setBaryonicFraction(double fb);
void setStellarRatio(double ratio);
void setStellarTotalMass(void);
void setGasDensity(void);
void setGasPower(double n);
