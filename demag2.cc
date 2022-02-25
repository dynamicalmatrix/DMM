#include <math.h>
#include "demag.h"


double Oc_Atan2(double y,double x)
{
#define Oc_Atan2(y,x) atan2((y),(x))
  if(x==0.0 && y==0.0) return 0.0;
  return atan2(y,x);
#undef Oc_Atan2
}

extern "C" {
static int AS_Compare(const void* px,const void* py)
{
  REALWIDE x=fabs(*((const REALWIDE *)px));
  REALWIDE y=fabs(*((const REALWIDE *)py));
  if(x<y) return 1;
  if(x>y) return -1;
  return 0;
}
}

static REALWIDE
AccurateSum(int n,REALWIDE *arr)
{
  qsort(arr,n,sizeof(REALWIDE),AS_Compare);

  REALWIDE sum,corr,y,u,t,v,z,x;
  sum=arr[0]; corr=0;
  for(int i=1;i<n;i++) {
    x=arr[i];
    y=corr+x;
    u=x-(y-corr);
    t=y+sum;
    v=y-(t-sum);
    z=u+v;
    sum=t+z;
    corr=z-(sum-t);
  }
  return sum;
}

static REALWIDE
SelfDemagDx(REALWIDE x,REALWIDE y,REALWIDE z)
{ 

  if(x==y && y==z) return 4*PI*x*x*x/3.;  

  REALWIDE xsq=x*x,ysq=y*y,zsq=z*z;
  REALWIDE diag=sqrt(xsq+ysq+zsq);
  REALWIDE arr[15];

  REALWIDE mpxy = (x-y)*(x+y);
  REALWIDE mpxz = (x-z)*(x+z);

  arr[0] = -4*(2*xsq*x-ysq*y-zsq*z);
  arr[1] =  4*(xsq+mpxy)*sqrt(xsq+ysq);
  arr[2] =  4*(xsq+mpxz)*sqrt(xsq+zsq);
  arr[3] = -4*(ysq+zsq)*sqrt(ysq+zsq);
  arr[4] = -4*diag*(mpxy+mpxz);

  arr[5] = 24*x*y*z*atan(y*z/(x*diag));
  arr[6] = 12*(z+y)*xsq*log(x);

  arr[7] = 12*z*ysq*log((sqrt(ysq+zsq)+z)/y);
  arr[8] = -12*z*xsq*log(sqrt(xsq+zsq)+z);
  arr[9] = 12*z*mpxy*log(diag+z);
  arr[10] = -6*z*mpxy*log(xsq+ysq);

  arr[11] =  12*y*zsq*log((sqrt(ysq+zsq)+y)/z);
  arr[12] = -12*y*xsq*log(sqrt(xsq+ysq)+y);
  arr[13] =  12*y*mpxz*log(diag+y);
  arr[14] =  -6*y*mpxz*log(xsq+zsq);

  REALWIDE Dx = AccurateSum(15,arr)/3.;

  return Dx;
}

static REALWIDE
f(REALWIDE x,REALWIDE y,REALWIDE z)
{ 

  x=fabs(x); REALWIDE xsq=x*x;
  y=fabs(y); REALWIDE ysq=y*y;
  z=fabs(z); REALWIDE zsq=z*z; 

  REALWIDE R=xsq+ysq+zsq;
  if(R<=0.0) return 0.0;
  else       R=sqrt(R);

  REALWIDE piece[8];
  int piececount=0;
  if(z>0.) { 
    REALWIDE temp1,temp2,temp3;
    piece[piececount++] = 2*(2*xsq-ysq-zsq)*R;
    if((temp1=x*y*z)>0.)
      piece[piececount++] = -12*temp1*atan2(y*z,x*R);
    if(y>0. && (temp2=xsq+zsq)>0.) {
      REALWIDE dummy = log(((y+R)*(y+R))/temp2);
      piece[piececount++] = 3*y*zsq*dummy;
      piece[piececount++] = -3*y*xsq*dummy;
    }
    if((temp3=xsq+ysq)>0.) {
      REALWIDE dummy = log(((z+R)*(z+R))/temp3);
      piece[piececount++] = 3*z*ysq*dummy;
      piece[piececount++] = -3*z*xsq*dummy;
    }
  } else {
    if(x==y) {
      const double K = -2.45981439737106805379;
      piece[piececount++] = K*xsq*x;
    } else {
      piece[piececount++] = 2*(2*xsq-ysq)*R;
      if(y>0. && x>0.)
	piece[piececount++] = -6*y*xsq*log((y+R)/x);
    }
  }

  return AccurateSum(piececount,piece)/12.;
}

extern "C"{
REALWIDE GetSDA00(REALWIDE x,REALWIDE y,REALWIDE z,
	 REALWIDE dx,REALWIDE dy,REALWIDE dz)
{ 
  REALWIDE result=0.;
  if(x==0. && y==0. && z==0.) {
    result = SelfDemagDx(dx,dy,dz);
  } else {
    REALWIDE arr[27];
    arr[ 0] = -1*f(x+dx,y+dy,z+dz);
    arr[ 1] = -1*f(x+dx,y-dy,z+dz);
    arr[ 2] = -1*f(x+dx,y-dy,z-dz);
    arr[ 3] = -1*f(x+dx,y+dy,z-dz);
    arr[ 4] = -1*f(x-dx,y+dy,z-dz);
    arr[ 5] = -1*f(x-dx,y+dy,z+dz);
    arr[ 6] = -1*f(x-dx,y-dy,z+dz);
    arr[ 7] = -1*f(x-dx,y-dy,z-dz);

    arr[ 8] =  2*f(x,y-dy,z-dz);
    arr[ 9] =  2*f(x,y-dy,z+dz);
    arr[10] =  2*f(x,y+dy,z+dz);
    arr[11] =  2*f(x,y+dy,z-dz);
    arr[12] =  2*f(x+dx,y+dy,z);
    arr[13] =  2*f(x+dx,y,z+dz);
    arr[14] =  2*f(x+dx,y,z-dz);
    arr[15] =  2*f(x+dx,y-dy,z);
    arr[16] =  2*f(x-dx,y-dy,z);
    arr[17] =  2*f(x-dx,y,z+dz);
    arr[18] =  2*f(x-dx,y,z-dz);
    arr[19] =  2*f(x-dx,y+dy,z);

    arr[20] = -4*f(x,y-dy,z);
    arr[21] = -4*f(x,y+dy,z);
    arr[22] = -4*f(x,y,z-dz);
    arr[23] = -4*f(x,y,z+dz);
    arr[24] = -4*f(x+dx,y,z);
    arr[25] = -4*f(x-dx,y,z);

    arr[26] =  8*f(x,y,z);

    result=AccurateSum(27,arr);
  }
  return result;
}
}

static REALWIDE
g(REALWIDE x,REALWIDE y,REALWIDE z)
{ 

  REALWIDE result_sign=1.0;
  if(x<0.0) result_sign *= -1.0;  if(y<0.0) result_sign *= -1.0;
  x=fabs(x); y=fabs(y); z=fabs(z);  

  REALWIDE xsq=x*x,ysq=y*y,zsq=z*z;
  REALWIDE R=xsq+ysq+zsq;
  if(R<=0.0) return 0.0;
  else       R=sqrt(R);

  REALWIDE piece[7];
  int piececount=0;
  piece[piececount++] = -2*x*y*R;;
  if(z>0.) { 
    piece[piececount++] = -z*zsq*atan2(x*y,z*R);
    piece[piececount++] = -3*z*ysq*atan2(x*z,y*R);
    piece[piececount++] = -3*z*xsq*atan2(y*z,x*R);

    REALWIDE temp1,temp2,temp3;
    if((temp1=xsq+ysq)>0.)
      piece[piececount++] = 6*x*y*z*log((z+R)/sqrt(temp1));

    if((temp2=ysq+zsq)>0.)
      piece[piececount++] = y*(3*zsq-ysq)*log((x+R)/sqrt(temp2));

    if((temp3=xsq+zsq)>0.)
      piece[piececount++] = x*(3*zsq-xsq)*log((y+R)/sqrt(temp3));

  } else {
    if(y>0.) piece[piececount++] = -y*ysq*log((x+R)/y);
    if(x>0.) piece[piececount++] = -x*xsq*log((y+R)/x);
  }

  return (1./6.)*result_sign*AccurateSum(piececount,piece);
}

extern "C"{
REALWIDE
GetSDA01(REALWIDE x,REALWIDE y,REALWIDE z,
	 REALWIDE l,REALWIDE h,REALWIDE e)
{ 

  REALWIDE arr[27];

  arr[ 0] = -1*g(x-l,y-h,z-e);
  arr[ 1] = -1*g(x-l,y-h,z+e);
  arr[ 2] = -1*g(x+l,y-h,z+e);
  arr[ 3] = -1*g(x+l,y-h,z-e);
  arr[ 4] = -1*g(x+l,y+h,z-e);
  arr[ 5] = -1*g(x+l,y+h,z+e);
  arr[ 6] = -1*g(x-l,y+h,z+e);
  arr[ 7] = -1*g(x-l,y+h,z-e);

  arr[ 8] =  2*g(x,y+h,z-e);
  arr[ 9] =  2*g(x,y+h,z+e);
  arr[10] =  2*g(x,y-h,z+e);
  arr[11] =  2*g(x,y-h,z-e);
  arr[12] =  2*g(x-l,y-h,z);
  arr[13] =  2*g(x-l,y+h,z);
  arr[14] =  2*g(x-l,y,z-e);
  arr[15] =  2*g(x-l,y,z+e);
  arr[16] =  2*g(x+l,y,z+e);
  arr[17] =  2*g(x+l,y,z-e);
  arr[18] =  2*g(x+l,y-h,z);
  arr[19] =  2*g(x+l,y+h,z);

  arr[20] = -4*g(x-l,y,z);
  arr[21] = -4*g(x+l,y,z);
  arr[22] = -4*g(x,y,z+e);
  arr[23] = -4*g(x,y,z-e);
  arr[24] = -4*g(x,y-h,z);
  arr[25] = -4*g(x,y+h,z);

  arr[26] =  8*g(x,y,z);

  return AccurateSum(27,arr);
}
}
