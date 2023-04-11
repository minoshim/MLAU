// 2D MHD simulation
// MPI+OpenMP

#include "global.hpp"
using namespace global;

#include "new_delete.hpp"
#include "init.hpp"
#include "dataio.hpp"
#include "cflcheck.hpp"
#include "bkup.hpp"

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

  new_delete();
  init_grid(mpi_rank);
  init_plasma(mpi_rank);

#if (CFLCHECK)
  cflcheck(&dt,!(mpi_rank));
#endif
  bkup_load(&n,&cnt,&tim,&dt,&trec,nx*ny,mpi_rank);
  if (n == 0){
    dataio(n,cnt,tim,mpi_rank);
  }
  
  stim=MPI_Wtime();
  while(n++ < nmax){
    tim+=dt;

    mhd_fd_ct_2d(ro,mx,my,mz,en,bx,by,bz,
		 dt,dx,dy,nx,ny,XOFF,YOFF,gam,
		 mpi_rank,mpi_numx,mpi_numy);

#if (CFLCHECK)
    cflcheck(&dt,0);
#endif

    // Output
    if ((n % nrec) == 0){  // Output at constant step
    // if (tim >= trec){	   // Output at constant period

      // NAN/INF check
      // MPI_Barrier(MPI_COMM_WORLD);
      // for (int ss=0;ss<nx*ny;ss++){
      // 	if (finite(ro[ss]) == 0){
      // 	  puts("Calculate halted.");
      // 	  MPI_Abort(MPI_COMM_WORLD,-1);
      // 	}
      // }

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
        new_delete();
        MPI_Finalize();
        return 0;
      }
    }

  }
  etim=MPI_Wtime();
  if (mpi_rank == 0) printf("%lu sec is required for computation.\n",(unsigned long )(etim-stim));

  new_delete();
  MPI_Finalize();
  return 0;
}
