/*  Given a hierarchical clustering through a sequence of agglomerations,
    derive the assignment for the top (lev-1)th level of the hierarchy. 
    
    Parameters:                                                    
    n:          number of observations                             
    ia,ib:      vectors of dimension N defining the agglomerations 
    lev:        number of clusters in partition.           
    iclass:     n-dimensional vector of cluster assignments       

    C code written by Ranjan Maitra, Baltimore (07/02/02) */ 

void hclass(int n, int *ia, int *ib, int lev, int *iclass)
{
  int i,j,k;
  for (i=0;i<n;i++) {
    iclass[i]=0;
  }
  j=lev-1;
  for (i=n-lev;i<(n-1);i++) {
    iclass[ib[i]]=j;
    for (k=(n-lev-1);k>=0;k--) {
      if (iclass[ia[k]]==j) {
	iclass[ib[k]]=j;
	}
    }
    j--;
  }
  j=0;
  iclass[ia[n-2]]=j;
  for (k=(n-lev-1);k>=0;k--) {
    if (iclass[ia[k]]==j) {
      iclass[ib[k]]=j;
    }
  }
  return;
}

