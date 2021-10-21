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
  time_t stim,etim;

  new_delete();
  init_grid();
  init_plasma();
  dataio(n,cnt,tim);

#if (CFLCHECK)
  cflcheck(dr,&dt,1);
#endif
  
  time(&stim);
  while(n++ < nmax){
    tim+=dt;

    setbc();
    mhd_fd_ct_3d(ro,mx,my,mz,en,bx,by,bz,
		 dt,dx,dy,dz,nx,ny,nz,XOFF,YOFF,ZOFF,gam);

#if (CFLCHECK)
    cflcheck(dr,&dt,0);
    // I confirm that dt does not change significantly till t=pi in OTvortex
#endif

    if ((n % nrec) == 0){
      cnt++;
      printf("%d / %d iterations finished.\n",n,nmax);
      dataio(n,cnt,tim);
    }
  }
  time(&etim);
  printf("%lu sec is required for computation.\n",etim-stim);

  new_delete();
  return 0;
}
