#include <stdio.h>
#include <math.h>

#define FROM -12.
#define TO 7.
#define STEPS 100

int main(void)
{
  int i;
  unsigned long int n_f, n_d;
  float number_f;
  double number_d, old, pos;

  for(i = 0; i < STEPS; i ++)
  {
    pos = pow(10, i * (TO - FROM) / STEPS + FROM);
    number_f = number_d = pos;
    n_f = *(unsigned int *) & number_f - 1;
    n_d = *(unsigned long int *) & number_d - 1;

    if (number_d != old)
    {
      old = number_d;
      printf("%.30f\t%.30f\t%.30f\n", pos, number_f - *(float *) & n_f, number_d - *(double *) & n_d);
    }
  }

  return 0;
}
