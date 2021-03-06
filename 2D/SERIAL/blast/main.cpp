#include "global.hpp"
using namespace global;

#include "new_vals.hpp"
#include "init_grid.hpp"
#include "init_plasma.hpp"
#include "dataio.hpp"
#include "setbc.hpp"
#include "delete_vals.hpp"
#include "cflcheck.hpp"
#include "cflcomment.hpp"

int main(void)
{
  int n=0,cnt=0;
  double tim=0;
  clock_t stim,etim;

  new_vals();
  init_grid();
  init_plasma();
  dataio(n,cnt,tim);

#if (CFLCHECK)
  cflcheck(dr,&dt);
#endif
  cflcomment(dr,dt);
  
  stim=clock();
  while(n++ < nmax){
    tim+=dt;

    setbc();
    mhd_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
    		 dt,dx,dy,nx,ny,XOFF,YOFF,gam);

#if (CFLCHECK)
    cflcheck(dr,&dt);
#endif
  
    for (int ss=0;ss<nx*ny;ss++){
      if (finite(ro[ss]) != 1){
	puts("Calculation halted");
	delete_vals();
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
  printf("%f sec is required for computation.\n",(double)(etim-stim)/CLOCKS_PER_SEC);

  delete_vals();
  return 0;
}
