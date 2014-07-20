#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cfacP.h"
#include "structure.h"

static char _ename[][3] = 
{"H", "He", "Li", "Be", "B", "C", "N", "O", "F",
 "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
 "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", 
 "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
 "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", 
 "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
 "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
 "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
 "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
 "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
 "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"};

static double _emass[] = 
{1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 
 40, 45, 48, 51, 52, 55, 56, 58, 59, 64, 65, 70, 73, 75, 79, 80, 84, 85,
 88, 89, 91, 93, 96, 98, 101, 103, 106, 108, 112, 115, 119, 122, 128, 127,
 131, 133, 137, 139, 140, 141, 144, 147, 150, 152, 157, 159, 162, 165, 167,
 169, 173, 175, 178, 181, 184, 186, 190, 190, 195, 197, 200, 204, 207, 209, 
 210, 210, 222, 223, 226, 227, 232, 231, 238, 237, 242, 243, 247, 247, 249, 
 254, 253, 256, 254, 257, 257, 260, 263, 262, 265, 266};


int cfac_set_atom(cfac_t *cfac, const char *s,
    double z, double mass, double rn) {
  cfac_nucleus_t *atom = &cfac->nucleus;
  unsigned int i, n_elements = sizeof(_emass)/sizeof(double);

  if (s == NULL) return -1;
  if (strlen(s) == 0) {
    if (z <= 0) {
      printf("atomic symbol and z cannot be both unset\n");
    }
    s = _ename[(int)(z-0.8)];
  }
  strncpy(atom->symbol, s, 2); 
  
  for (i = 0; i < n_elements; i++) {
    if (strncasecmp(_ename[i], s, 2) == 0) {
      if (z <= 0) atom->atomic_number = i+1;
      if (mass <= 0) atom->mass = _emass[i];
      break;
    }
  }
  if (i == n_elements) return -1;
  
  cfac->anum = i + 1;

  if (z > 0) {
    atom->atomic_number = z;
  } 
  if (mass > 0.0) {
    atom->mass = mass;
  }
  if (rn < 0.0) {
    atom->rn = 2.2677E-5 * pow(atom->mass, 1.0/3);
  } else {
    atom->rn = rn;
  }
  
  /* allocate & init per-charge-state level arrays */
  cfac->levels_per_ion = malloc(sizeof(ARRAY)*(cfac->anum + 1));

  for (i = 0; i <= cfac->anum; i++) {
    ArrayInit(&cfac->levels_per_ion[i], sizeof(LEVEL_ION), 512, NULL, NULL);
  }


  return 0;
}

double cfac_get_atomic_mass(const cfac_t *cfac) {
  return cfac->nucleus.mass;
}

double cfac_get_atomic_number(const cfac_t *cfac) {
  return cfac->nucleus.atomic_number;
}

const char *cfac_get_atomic_symbol(const cfac_t *cfac) {
  return cfac->nucleus.symbol;
}

double cfac_get_atomic_effective_z(const cfac_t *cfac, double r) {
  cfac_nucleus_t atom = cfac->nucleus;
  double x, y;
  if (r > atom.rn) {
    return (double) atom.atomic_number;
  } else {
    x = r/atom.rn;
    y = 3.0 - x*x;
    y = x*y*0.5*(atom.atomic_number);
    return y;
  }
}

double cfac_get_atomic_rn(const cfac_t *cfac) {
  return cfac->nucleus.rn;
}
