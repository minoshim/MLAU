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
  time_t stim,etim;

  new_vals();
  init_grid();
  init_plasma();
  dataio(n,cnt,tim);

#if (CFLCHECK)
  cflcheck(dr,&dt);
#endif
  cflcomment(dr,dt);
  
  time(&stim);
  while(n++ < nmax){
    tim+=dt;

    setbc();
    mhd_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
    		 dt,dx,dy,nx,ny,XOFF,YOFF,gam);

#if (CFLCHECK)
    cflcheck(dr,&dt);
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

  delete_vals();
  return 0;
}
