// 2D MHD simulation
// MPI+OpenMP

#include "global.hpp"
using namespace global;

#include "new_vals.hpp"
#include "init_grid.hpp"
#include "init_plasma.hpp"
#include "dataio.hpp"
#include "delete_vals.hpp"
#include "cflcheck.hpp"
#include "cflcomment.hpp"
#include "bkup_save.hpp"
#include "bkup_load.hpp"

int main(int argc, char* argv[])
{
  // MPI Send/Recv arguments
  int mpi_rank;
  int mpi_num;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
  if (mpi_num != MNP){
    if (mpi_rank == 0) puts("MPI number of process is not correct. Check MNP in global.hpp");
    MPI_Finalize();
    return 0;
  }

  int n=0,cnt=0;
  double tim=0;
  time_t stim,etim;
  double tlimit=1400;           // Calculation time limit in min

  new_vals();
  init_grid(mpi_rank);
  init_plasma(mpi_rank);

#if (CFLCHECK)
  cflcheck(&dt);
#endif
  bkup_load(&n,&cnt,&tim,&dt,&trec,nx*ny,mpi_rank);
  if (n == 0){
    dataio(n,cnt,tim,mpi_rank);
  }
  cflcomment(dr,dt,mpi_rank);
  
  stim=MPI_Wtime();
  while(n++ < nmax){
    tim+=dt;

    mhd_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
		 dt,dx,dy,nx,ny,XOFF,YOFF,gam,
		 mpi_rank,mpi_numx,mpi_numy);

#if (DIFF)
    if (eta0 != 0){
      diff_ctfield_e(en,bx,by,bz,eta0,
		     dt,dx,dy,nx,ny,XOFF,YOFF,gam,
		     mpi_rank,mpi_numx,mpi_numy);
    }
    if (nu0  != 0){
      ns_viscous_2d(ro,mx,my,mz,en,nu0,
		    dt,dx,dy,nx,ny,XOFF,YOFF,gam,
		    mpi_rank,mpi_numx,mpi_numy);
    }
#endif
    
#if (CFLCHECK)
    cflcheck(&dt);
#endif

    // Output
    // if ((n % nrec) == 0){  // Output at constant step
    if (tim >= trec){	   // Output at constant period

      // NAN/INF check
      MPI_Barrier(MPI_COMM_WORLD);
      for (int ss=0;ss<nx*ny;ss++){
	if (finite(ro[ss]) == 0){
	  puts("Calculate halted.");
	  MPI_Abort(MPI_COMM_WORLD,-1);
	}
      }

      trec+=dtrec;
      cnt++;
      bkup_save(n,cnt,tim,dt,trec,nx*ny,mpi_rank);
      dataio(n,cnt,tim,mpi_rank);
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0) printf("%d / %d iterations finished.\n",n,nmax);

      // Elapse time check
      MPI_Barrier(MPI_COMM_WORLD);
      if (((MPI_Wtime()-stim)/60.) > tlimit){
        if (mpi_rank == 0)
          printf("Exceed time limit %f min. The job is terminated.\n",tlimit);
        delete_vals();
        MPI_Finalize();
        return 0;
      }
    }

  }
  etim=MPI_Wtime();
  if (mpi_rank == 0) printf("%lu sec is required for computation.\n",(unsigned long )(etim-stim));

  delete_vals();
  MPI_Finalize();
  return 0;
}
