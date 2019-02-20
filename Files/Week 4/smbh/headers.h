#define PI acos(0) * 2
#define LN_LAMBDA 2.3

#define h 0.678
#define G0 0.449849
// #define G0 4.49849e-06 // kpc gyr m0
#define H 0.069158

void setR_vir(double r);
void printConstants(void);
double gasDensity(double r);
double getNorm(double *vector);
double darkMatterMass(double r);
double darkMatterDensity(double r);
double darkMatterDensity0(double c);
double getLocalSoundSpeed(double z);
double gravitationalForce(double r);
double baryonicMassHernquist(double r);
double baryonicDensityHernquist(double r);
void baseCase(struct reb_simulation* sim);
double SMBHAccretion(double *position, double *speed);
double *dynamicalFrictionDM(double *position, double *speed);
double *dynamicalFrictionGas(double *position, double *speed);

void sphericalToCartesian(double *r, double *theta, double *phi);
void printStatus(struct reb_simulation *sim, const char *filename);

struct reb_simulation* setupSimulation(double mass, double *position, double *speed,
                            void (*additional_force)(struct reb_simulation*),
                            int velocity_dependent);
void runSimulation(struct reb_simulation *sim, int n_points, int save_points, const char *filename);

void run(double *positions, double *speeds, double smbh_mass, double dt, int n_points, int n_saves, const char *filename);
