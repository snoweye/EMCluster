#define inf 1e+40;
#include<stdlib.h>
#include<math.h>
#include "array.h"

/*  HIERARCHICAL CLUSTERING using (user-specified) criterion. 

    Parameters:                                               
    n,m               number of observations and coordinates  
    data(n,m)         input data matrix,                      
    iopt              clustering criterion to be used,        
    IA, IB, CRIT      history of agglomerations; dimension N-1 each       

    C code rewritten by Ranjan Maitra, Baltimore (07/02/02) from f77
    code of F. Murtagh, ESA/ESO/STECF, Garching (1986). Note that I use the 
    true distances, not squared Euclidean which FM used. */


int ioffset(int n, int i, int j)
{                               /*  Map (i,j)th element of upper half */
  return(j+i*n-(i+1)*(i+2)/2);  /*  diagonal symmetric matrix to vector */
} 

int imin(int aa, int bb)
{
  if (aa<bb) {
    return(aa);
  }
  else return(bb);
}

int imax(int aa, int bb)
{
  if (aa<bb) {
    return(bb);
  }
  else return(aa);
}

double fmin(double aa, double bb)
{
  if (aa<bb) {
    return(aa);
  }
  else return(bb);
}

double fmax(double aa, double bb)
{
  if (aa<bb) {
    return(bb);
  }
  else return(aa);
}

void hc(int n, int m, int iopt, double **data, int *ia, int *ib, double *crit)
{
  int i,i2,ind,ind1,ind2,ind3,ncl,j,j2,k,len,iflag,*nn,*flag;
  double dmin,r1,r2,x,xx,*disnn,*membr,*diss,*disso;

  ncl = n;
  len=n*(n-1)/2;
  
  MAKE_VECTOR(flag,n);
  MAKE_VECTOR(nn,n);
  MAKE_VECTOR(disnn,n);
  MAKE_VECTOR(membr,n);
  MAKE_VECTOR(disso,len);
  MAKE_VECTOR(diss,len);

  for (i=0; i<n; ++i) {
    nn[i]=0;
    disnn[i]=0;
    membr[i] = 1;
    flag[i] = 1;
  }

  for (i=0; i<len; ++i) {
    disso[i]=0;
    diss[i]=0;
  }

  for (i=0;i<(n-1);i++) {   /*  Construct dissimilarity matrix */
    for (j=(i+1);j<n;j++) {
      ind = ioffset(n, i, j);
      for (k=0;k<m;k++) { /* Computing Euclidean distance */
	diss[ind]+=(data[i][k]-data[j][k])*(data[i][k]-data[j][k]);
      }
      if (iopt == 1) {  /*for the case of the min. var. method where merging */
	diss[ind] /= 2.; /*criteria are defined in terms of variances */ 
      }
      else {
	diss[ind]=sqrt(diss[ind]); /* true Euclidean distance */
      }
      disso[ind] = diss[ind];
    }
  }
  
  for (i=0;i<(n-1);++i) { /*Carry out an agglomeration - first create list of NNs */
    int jm;
    dmin = inf;
    jm=i+1;
    for (j=(i+1);j<n;++j) {
      ind=ioffset(n,i,j);
      if (diss[ind]<dmin) {
	dmin = diss[ind];
	jm = j;
      }
    }
    nn[i]=jm;
    disnn[i]=dmin;
  }

  while (ncl>1) { /*     Next, determine least diss. using list of NNs */ 
    int jj=0,jm=0,im=0;
    dmin=inf;
    for (i= 0;i<(n-1);++i) {
      if (flag[i]==1) {
	if (disnn[i]<dmin) {
	  dmin=disnn[i];
	  im=i;
	  jm=nn[i];
	  }
      }
    }

    i2=imin(im,jm); /*  This allows an agglomeration to be carried out. */
    j2=imax(im,jm);
    ia[n-ncl]=i2;
    ib[n-ncl]=j2;
    crit[n-ncl]=dmin;
    
    flag[j2]=0;     /*  Update dissimilarities from new cluster. */
    dmin=inf;
    for (k=0;k<n;k++) {
      if ((flag[k]==1) && (k!=i2)) {
	x=membr[i2]+membr[j2]+membr[k];
	ind1=ioffset(n,imin(i2,k),imax(i2,k));
	ind2=ioffset(n,imin(j2,k),imax(j2,k));
	ind3=ioffset(n,i2,j2);
	xx=diss[ind3];
	
	if (iopt == 1) { /*  WARD'S MINIMUM VARIANCE METHOD - IOPT=1. */
	  diss[ind1] = (membr[i2] + membr[k]) * diss[ind1] + 
	    (membr[j2] + membr[k]) * diss[ind2] - membr[k] * xx;
	  diss[ind1] /= x;
	}
	
	if (iopt == 2) { /*  SINGLE LINK METHOD - IOPT=2. */ 
	  r1 = diss[ind1]; /* Computing MIN */
	  r2 = diss[ind2];
	  diss[ind1]=fmin(r1,r2);
	}
	
	if (iopt==3) {  /*  COMPLETE LINK METHOD - IOPT=3. */
	  r1=diss[ind1]; /* Computing MAX */
	  r2=diss[ind2];
	  diss[ind1] = fmax(r1,r2);
	}
	
	if (iopt==4) { /* AVG. LINK (OR GROUP AVERAGE) METHOD - IOPT=4. */
	  diss[ind1] = (membr[i2]*diss[ind1]+membr[j2]*diss[ind2])/
	    (membr[i2]+membr[j2]); 
	}
	
	if (iopt == 5) { /*  MCQUITTY'S METHOD - IOPT=5. */
	    diss[ind1] = diss[ind1]*0.5 + diss[ind2]*0.5;
	}
	
	if (iopt==6) { /*  MEDIAN (GOWER'S) METHOD - IOPT=6. */
	  diss[ind1]=diss[ind1]*0.5+diss[ind2]*0.5-xx*0.25;
	}
	
	if (iopt==7) { /*  CENTROID METHOD - IOPT=7. */
	  diss[ind1] = (membr[i2]*diss[ind1]+membr[j2]*diss[ind2] - 
			membr[i2]*membr[j2]*xx/(membr[i2]+membr[j2])
			)/(membr[i2] + membr[j2]);
	}
	if (i2<=k) {
	  if (diss[ind1]<dmin) { 
	    dmin=diss[ind1];
	    jj=k;
	  }
	}
      }
    }
    membr[i2]+=membr[j2];
    disnn[i2]=dmin;
    nn[i2] = jj;
    for (i=0;i<(n-1);++i) { /* Update list of NNs insofar as  required. */
      iflag=0;
      if (iflag==0) {
	if (flag[i]==0) {
	  iflag=1;
	}
	if (iflag==0) {
	  if ((nn[i]!=i2) && (nn[i]!=j2))  { 
	    iflag=1;
	  }
	  if (iflag==0) { /*  (Redetermine NN of i) */
	    int jj=0;
	    dmin = inf;
	    for (j=(i+1);j<n;++j) {
	      ind = ioffset(n,i,j);
	      if ((flag[j]==1) & (diss[ind]<=dmin)) {
		dmin = diss[ind];
		jj = j;
	      }
	    }
	    nn[i] = jj;
	    disnn[i] = dmin;
	  }
	}
      }
    }
    ncl--; /*  Repeat previous steps until n-1 agglomerations carried out. */
  }
  FREE_VECTOR(flag);
  FREE_VECTOR(nn);
  FREE_VECTOR(membr);
  FREE_VECTOR(disnn);
  FREE_VECTOR(disso);
  FREE_VECTOR(diss);
  return;
}
