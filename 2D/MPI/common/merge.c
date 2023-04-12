#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
/* Merge 2D MPI-decomposed data */
/* > ./merge.out FILE_DIRECTORY OUTPUT_DIRECTORY */
{
  if (argc !=3){
    fputs("Set file directory and output directory\n",stderr);
    fputs("Command Syntax: ./a.out FILDIR OUTDIR\n",stderr);
    return 0;
  }
  char *chtmp;
  char fildir[100];
  char outdir[100];
  chtmp=argv[1];
  strcpy(fildir,chtmp);
  strcat(fildir,"/");
  chtmp=argv[2];
  strcpy(outdir,chtmp);
  strcat(outdir,"/");
  printf("File directory is %s\n",fildir);
  printf("Output directory is %s\n",outdir);
  puts("");

  const int nstt=0;		/* Start time index */
  const int smax=8;		/* Number of MHD variables */
  int i,j,m,n,mx,my,s,cnt;
  int nx=0,ny=0,nt=0;
  int nxall=0,nyall=0;
  int xoff=0,yoff=0;
  int mpix=0,mpiy=0;
  int itmp;

  FILE *infil,*outfil;
  char filname[100];
  double dummy;

  /* Read offset data */
  sprintf(filname,"%soffsets.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    int tmp[2]={0};
    itmp=fscanf(infil,"%d %d\n",&tmp[0],&tmp[1]);
    xoff=tmp[0];
    yoff=tmp[1];
    fclose(infil);

    puts("Offset data is found.");
    printf("Offsets in X and Y directions are %d and %d.\n",xoff,yoff);
    puts("");
  } else{
    puts("No offset data is found.");
    return 0;
  }
  
  /* Check X-array data */
  cnt=0;
  sprintf(filname,"%sx_%05d.dat",fildir,cnt++);
  if ((infil=fopen(filname,"r"))!=NULL){

    while(fscanf(infil,"%lf",&dummy) !=EOF) nx++;
    mpix++;
    nxall+=(nx-2*xoff);
    fclose(infil);

    sprintf(filname,"%sx_%05d.dat",fildir,cnt++);
    while((infil=fopen(filname,"r"))!=NULL){
      mpix++;
      nxall+=(nx-2*xoff);
      sprintf(filname,"%sx_%05d.dat",fildir,cnt++);
      fclose(infil);
    }

    puts("X-array data is found.");
    printf("X grid number per 1 MPI process is %d.\n",nx);
    printf("MPI process number along X is %d.\n",mpix);
    printf("Total X grid (w.o. offset) is %d\n",nxall);
    puts("");
  } else{
    puts("No X-array data is found.");
    return 0;
  }
  /* Read, merge, and write X-array data (removing offset) */
  double x[nxall];
  for (m=0;m<mpix;m++){
    double xtmp[nx];
    sprintf(filname,"%sx_%05d.dat",fildir,m);
    infil=fopen(filname,"r");
    i=0;
    while(fscanf(infil,"%lf",&dummy) != EOF) xtmp[i++]=dummy;
    fclose(infil);
    for (i=0;i<nx-2*xoff;i++){
      x[m*(nx-2*xoff)+i]=xtmp[i+xoff];
    }
  }
  sprintf(filname,"%smerge_x.dat",outdir);
  outfil=fopen(filname,"w");
  for (i=0;i<nxall;i++) fprintf(outfil,"%f\n",x[i]);
  fclose(outfil);

  /* Check Y-array data */
  cnt=0;
  sprintf(filname,"%sy_%05d.dat",fildir,cnt++);
  if ((infil=fopen(filname,"r"))!=NULL){

    while(fscanf(infil,"%lf",&dummy) !=EOF) ny++;
    mpiy++;
    nyall+=(ny-2*yoff);
    fclose(infil);

    sprintf(filname,"%sy_%05d.dat",fildir,cnt++);
    while((infil=fopen(filname,"r"))!=NULL){
      mpiy++;
      nyall+=(ny-2*yoff);
      sprintf(filname,"%sy_%05d.dat",fildir,cnt++);
      fclose(infil);
    }

    puts("Y-array data is found.");
    printf("Y grid number per 1 MPI process is %d.\n",ny);
    printf("MPI process number along Y is %d.\n",mpiy);
    printf("Total Y grid (w.o. offset) is %d\n",nyall);
    puts("");
  } else{
    puts("No Y-array data is found.");
    return 0;
  }
  /* Read, merge, and write Y-array data (removing offset) */
  double y[nyall];
  for (m=0;m<mpiy;m++){
    double ytmp[ny];
    sprintf(filname,"%sy_%05d.dat",fildir,m);
    infil=fopen(filname,"r");
    j=0;
    while(fscanf(infil,"%lf",&dummy) != EOF) ytmp[j++]=dummy;
    fclose(infil);
    for (j=0;j<ny-2*yoff;j++){
      y[m*(ny-2*yoff)+j]=ytmp[j+yoff];
    }
  }
  sprintf(filname,"%smerge_y.dat",outdir);
  outfil=fopen(filname,"w");
  for (j=0;j<nyall;j++) fprintf(outfil,"%f\n",y[j]);
  fclose(outfil);

  /* Check Time-array data */
  sprintf(filname,"%st.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    while(fscanf(infil,"%lf",&dummy) != EOF) nt++;
    fclose(infil);

    puts("Time-array data is found");
    printf("Number of time step is %d.\n",nt);
    puts("");
  } else{
    puts("No Time-array data is found.");
    return 0;
  }

  /* Copy fildir/t.dat to outdir */
  char command[100];
  sprintf(command,"cp %st.dat %s",fildir,outdir);
  itmp=system(command);
  /* Copy fildir/params.dat to outdir */
  sprintf(command,"cp %sparams.dat %s",fildir,outdir);
  itmp=system(command);
  
  /* Data read, merge, write */
  int mpi=mpix*mpiy,flag;
  float *val[smax],*merge_val[smax];
  for (s=0;s<smax;s++){
    val[s]=(float*)malloc(sizeof(float)*nx*ny);
    merge_val[s]=(float*)malloc(sizeof(float)*nxall*nyall);
  }
  
  for (n=nstt;n<nt;n++){
    
    printf("Data read, merge, write at t = %d...\n",n);

    /* Initialize */
    for (s=0;s<smax;s++){
      for (j=0;j<nyall;j++){
	for (i=0;i<nxall;i++) merge_val[s][nxall*j+i]=0;
      }
    }
    /* Read and merge */
    flag=0;    
    for (m=0;m<mpi;m++){
      mx=m%mpix;
      my=m/mpix;

      sprintf(filname,"%soutdat_%05d_%05d.dat",fildir,n,m);
      if ((infil=fopen(filname,"rb"))!=NULL){
	flag++;
	/* Read */
	for (s=0;s<smax;s++){
	  itmp=fread(val[s],sizeof(float),nx*ny,infil);
	  for (j=0;j<ny-2*yoff;j++){
	    for (i=0;i<nx-2*xoff;i++){
	      merge_val[s][nxall*(my*nyall/mpiy+j)+(mx*nxall/mpix+i)]=val[s][nx*(j+yoff)+(i+xoff)];
	    }
	  }
	}
	fclose(infil);
      }
    }
    /* Write */
    if (flag == mpi){
      sprintf(filname,"%smerge_outdat_%05d.dat",outdir,n);
      outfil=fopen(filname,"wb");
      for (s=0;s<smax;s++) fwrite(merge_val[s],sizeof(float),nxall*nyall,outfil);
      fclose(outfil);
    }
  }

  for (s=0;s<smax;s++){
    free(val[s]);
    free(merge_val[s]);
  }
  return 0;
}
