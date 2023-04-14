#include "global.hpp"
using namespace global;

#include "new_delete.hpp"
#include "init.hpp"
#include "dataio.hpp"
#include "setbc.hpp"
#include "cflcheck.hpp"

int main(void)
{
  int n=0,cnt=0;
  double tim=0;
  clock_t stim,etim;

  new_delete();
  init_grid();
  init_plasma();
  dataio(n,cnt,tim);

#if (CFLCHECK)
  cflcheck(dr,&dt,1);
#endif
  
  stim=clock();
  while(n++ < nmax){
    tim+=dt;

    setbc();
    mhd_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
    		 dt,dx,dy,nx,ny,XOFF,YOFF,gam);

#if (CFLCHECK)
    cflcheck(dr,&dt,0);
#endif
  
    for (int ss=0;ss<nx*ny;ss++){
      if (finite(ro[ss]) != 1){
	puts("Calculation halted");
	new_delete();
	return 0;
      }
    }

    if (tim >= trec){
      cnt++;
      trec+=dtrec;
      printf("%d / %d iterations finished.\n",n,nmax);
      dataio(n,cnt,tim);
    }
  }
  etim=clock();
  printf("CPU time = %f sec.\n",(double)(etim-stim)/CLOCKS_PER_SEC);

  new_delete();
  return 0;
}
