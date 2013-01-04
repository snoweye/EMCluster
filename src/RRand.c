void RRand(int *N, int *TRUK,int *PREDK,int *trcl, int *prcl,
	   double *Rand,double *adjRand,double *Eindex)
{
  int i,j,n[(*TRUK)][(*PREDK)];
  double sumtr[(*TRUK)],sumpr[(*PREDK)],sumprsq,sumtrsq,sumsq,discordant,
    sumtrprsq;
  double term1, term2, term3;
  double nij2sum, nidot2sum, ndotj2sum, Wallace;

  for (i=0;i<(*TRUK);i++)    {
    for (j=0;j<(*PREDK);j++)     {
      n[i][j]=0;
    }
  }
  
  for (i=0;i<(*N);i++) {
    n[trcl[i]][prcl[i]]+=1;
  }

  sumtrsq=0.;
  for (i=0;i<(*TRUK);i++) {
    sumtr[i]=0.;
    for (j=0;j<(*PREDK);j++) {
      sumtr[i]+=n[i][j];    }
    sumtrsq+=sumtr[i]*sumtr[i];
  }
  
  sumprsq=0.;
  for (j=0;j<(*PREDK);j++) {
    sumpr[j]=0.;
    for (i=0;i<(*TRUK);i++) {
      sumpr[j]+=(double)n[i][j];    }
    sumprsq+=sumpr[j]*sumpr[j];
  }

  sumtrprsq=0.;
  for (i=0;i<(*TRUK);i++) {
    for (j=0;j<(*PREDK);j++) {
      sumtrprsq+=sumtr[i]*sumtr[i]*sumpr[j]*sumpr[j];
    }
  }

  (*Eindex)=sumtrprsq/((*N)*((double)(*N)-1) + (*N)*(double)(*N)/((*N)-1)) - (sumprsq + sumtrsq)/((*N)-1);
  (*Eindex)*=2.;
  (*Eindex)/=(*N)*((double)(*N)-1);
  
  sumsq=0.;
  for (i=0;i<(*TRUK);i++)    {
    for (j=0;j<(*PREDK);j++) {
      sumsq+=(double)n[i][j]*n[i][j];
    }
  }

  nij2sum=0.;
  for (i=0;i<(*TRUK);i++) {
    for (j=0;j<(*PREDK);j++) {
      nij2sum+=(double)n[i][j]*(n[i][j]-1)/2.0;
    }
  }

  nidot2sum=0.;
  for (i=0;i<(*TRUK);i++) {
    nidot2sum+=(double)sumtr[i]*(sumtr[i]-1)/2.0;
  }

  ndotj2sum=0.;
  for (i=0;i<(*PREDK);i++) {
    ndotj2sum+=(double)sumpr[i]*(sumpr[i]-1)/2.0;
  }

  Wallace=nij2sum/nidot2sum;
  discordant = 0.5*(sumtrsq + sumprsq) - sumsq ;

  (*Rand)=1.0-discordant/((double)(*N)*((double)(*N)-1.)/2.);

  term3 = nidot2sum * ndotj2sum / ((double)(*N)*((double)(*N)-1.)/2.);

  term1 = nij2sum - term3;

  term2 = (nidot2sum + ndotj2sum)/2 - term3;

  (*adjRand)= term1/term2;

}

