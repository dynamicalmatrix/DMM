/* FORTRAN-C++ interface */

#include "demag.h"

//float tdemag(int *ex, int *ey, int *ez, int *erow,
float tdemag_(int *ex, int *ey, int *ez, int *erow,
             int *ecolumn, float *Cellsize, float *ZCellsize)
{
  float cs, zcs;
  REALWIDE x, y, z, dx, dy, dz;
  REALWIDE GetSDA00(REALWIDE,REALWIDE,REALWIDE,REALWIDE,REALWIDE,REALWIDE);
  REALWIDE GetSDA01(REALWIDE,REALWIDE,REALWIDE,REALWIDE,REALWIDE,REALWIDE);

  cs = *Cellsize;
  zcs = *ZCellsize;

  // Unit conversion (from cm to m)
  dx=cs/100;
  dy=dx;
  dz=zcs/100;
  x = *ex*dx;
  y = *ey*dy;
  z = *ez*dz;

  switch (*erow) {
  case 1:
   // Nxx:
   if (*ecolumn==1) return (float)GetSDA00(x,y,z,dx,dy,dz)/dx/dy/dz;
   //Nxy
   if (*ecolumn==2) return (float)GetSDA01(x,y,z,dx,dy,dz)/dx/dy/dz;
   //Nxz
   if (*ecolumn==3) return (float)GetSDA01(x,z,y,dx,dz,dy)/dx/dy/dz;
   break;
  case 2:
   // Nyx:
   if (*ecolumn==1) return (float)GetSDA01(x,y,z,dx,dy,dz)/dx/dy/dz;
   // Nyy:
   if (*ecolumn==2) return (float)GetSDA00(y,x,z,dy,dx,dz)/dx/dy/dz;
   // Nyz:
   if (*ecolumn==3) return (float)GetSDA01(y,z,x,dy,dz,dx)/dx/dy/dz;
   break;
  case 3:
   // Nzx:
   if (*ecolumn==1) return (float)GetSDA01(x,z,y,dx,dz,dy)/dx/dy/dz;
   // Nzy:
   if (*ecolumn==2) return (float)GetSDA01(y,z,x,dy,dz,dx)/dx/dy/dz;
   //Nzz
   if (*ecolumn==3) return (float)GetSDA00(z,y,x,dz,dy,dx)/dx/dy/dz;
  }
  fprintf (stderr, "demag: wrong row/column!\n");
  exit (-1);
}
