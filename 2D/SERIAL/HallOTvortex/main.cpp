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
  int nw,nwmax;

  new_vals();
  init_grid();
  init_plasma();
  setbc();
  dataio(n,cnt,tim);

#if (CFLCHECK)
  cflcheck(dr,&dt, 0);		// MHD part
  cflcheck(dr,&dtw,1);		// Hall part
#endif
  cflcomment(dr,dt, 0);
  cflcomment(dr,dtw,1);

  time(&stim);
  while((n++ < nmax) && (tim < tmax)){
    tim+=dt;

    // Sub-cycling the Hall term
#if (HALL)
    nwmax=1+(int)(dt/dtw);
    for (nw=0;nw<nwmax;nw++){
      setbc();
      hall_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
		    dt/nwmax,dx,dy,nx,ny,XOFF,YOFF,eta_h);
    }
#endif
    
    setbc();
    mhd_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
    		 dt,dx,dy,nx,ny,XOFF,YOFF,gam);

#if (CFLCHECK)
    cflcheck(dr,&dt, 0);	// MHD part
    cflcheck(dr,&dtw,1);	// Hall part
#endif

    if (tim >= trec){
      trec+=dtrec;
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
