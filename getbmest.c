#include <R.h>
#include <Rmath.h>
#include <stdio.h>

void getbmest(double *no3lst, double *pNi, double *abmi, int *plabmi, double *abm_esti) {
  int labmi=*plabmi;
  int ord[labmi];
  int subsi[labmi];
  int i, j;
  double sumpni;
  int tmp;

  //Get min R*
  double minno3lst=no3lst[0];
  for(i=1; i<labmi; i++) {
    if((no3lst[i]<minno3lst)) {
      minno3lst=no3lst[i];
    }
  }

  //Set abm for min R*
  for(i=0; i<labmi; i++) {
    if(no3lst[i]==minno3lst) {
      subsi[i] = 1;
      abm_esti[i]=abmi[i];
    }
  }
  
  //nitrogen competition
  for(i=0; i<labmi; i++) {
    ord[i] = i;
  }
  for(i=0; i<labmi; i++) {
    for(j=i+1; j<labmi; j++) {
      if(no3lst[ord[i]] > no3lst[ord[j]]) {
        tmp=ord[i];
        ord[i]=ord[j];
        ord[j]=tmp;
      }
    }
  }

  for(i=1; i<labmi; i++) {
    sumpni=0;
    for(j=0; j<i; j++) {
      sumpni=sumpni+(pNi[ord[j]]*abm_esti[ord[j]])/pNi[ord[i]];
    }
    abm_esti[ord[j]]=abmi[ord[j]]-sumpni;
    if(abm_esti[ord[j]]<0) {
      abm_esti[ord[j]]=0;
    }
  }
}



