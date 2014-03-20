#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Print.h>
#include <Rmath.h>
#include <math.h>

// new functions based on the EM algorithm
SEXP ScanIGSGridCumSumNewC(SEXP ATS, SEXP gridCurS) {
  double *AT = REAL(ATS);
  double *gridCur = REAL(gridCurS);
  long long gridCurLen = length(gridCurS);
  SEXP ATCumSum;
  PROTECT(ATCumSum = allocVector(REALSXP, gridCurLen-1));
  double *ATCumSumPtr = REAL(ATCumSum);
  long long i, j;
  for (i=0; i<gridCurLen-1; i++){
    ATCumSumPtr[i] = 0;
    for (j=gridCur[i]-1; j<gridCur[i+1]-1; j++){
      ATCumSumPtr[i] += AT[j];
    }
  }
  UNPROTECT(1);
  return(ATCumSum);
}

SEXP GetP(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP errorS, SEXP maxIterS, SEXP pOriS) {
  double *AT = REAL(ATS);
  double *BT = REAL(BTS); 
  double *AN = REAL(ANS);
  double *BN = REAL(BNS);
  double *pOri = REAL(pOriS);
  double error = REAL(errorS)[0];
  double maxIter = REAL(maxIterS)[0];
  long long len = length(ATS);
  SEXP rS, pS;
  PROTECT(pS = allocVector(REALSXP, 2));
  PROTECT(rS = allocVector(REALSXP, len));
  
  double *r = REAL(rS);
  double *p = REAL(pS);
  double nIter = 0;
  double curError = 1;
  double pa, pb, paOld, pbOld, pa1, pa2, pb1, pb2;
  pa = pOri[0];
  pb = pOri[1];
  while (curError>error && nIter<maxIter){
    paOld = pa;
    pbOld = pb;
    for (long long i=0; i<len; i++){
      double temp = (AT[i]-BT[i])*log(pb/pa) + (AN[i]-BN[i])*log((1-pb)/(1-pa));
      if (temp>100){
        r[i] = exp(-temp);
      }else{
        r[i] = 1.0/(1.0+exp(temp));
      }
        // r[i] = 1.0/(1.0+pow(pb/pa, AT[i]-BT[i])*pow((1-pb)/(1-pa), AN[i]-BN[i]));
    }
    pa1 = pa2 = pb1 = pb2 = 0.0;
    for (long long i=0; i<len; i++){
      pa1 += AT[i]*r[i] + BT[i]*(1-r[i]);
      pa2 += (AT[i]+AN[i])*r[i] + (BT[i]+BN[i])*(1-r[i]);
      pb1 += AT[i]*(1-r[i]) + BT[i]*r[i];
      pb2 += (AT[i]+AN[i])*(1-r[i]) + (BT[i]+BN[i])*r[i];
    }
    pa = pa1/pa2;
    pb = pb1/pb2;
    curError = sqrt(pow(pa-paOld,2) + pow(pb-pbOld,2));
    nIter = nIter + 1;
  }
  p[0] = pa;
  p[1] = pb;
  
  UNPROTECT(2);
  return(pS);
}

SEXP Lik(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP pS) {
  double *AT = REAL(ATS);
  double *BT = REAL(BTS);
  double *AN = REAL(ANS);
  double *BN = REAL(BNS);
  double *p = REAL(pS);
  double pa = p[0];
  double pb = p[1];
  long long len = length(ATS);
  SEXP likS;
  PROTECT(likS = allocVector(REALSXP, 1));
  double *lik = REAL(likS);
  double likVal = 0;
  for (long long i=0; i<len; i++){
    // likVal += log(pow(pa,AT[i])*pow(1-pa,AN[i])*pow(pb,BT[i])*pow(1-pb,BN[i]) + pow(pa,BT[i])*pow(1-pa,BN[i])*pow(pb,AT[i])*pow(1-pb,AN[i])); // has problem wen AT too large
    if (pa*(1-pa)*pb*(1-pb) !=0){
      double temp = (BT[i]-AT[i])*log(pa/pb) + (BN[i]-AN[i])*log((1-pa)/(1-pb));
      if (temp<100){
        likVal += AT[i]*log(pa) + AN[i]*log(1-pa) + BT[i]*log(pb) + BN[i]*log(1-pb) + log(1 + exp(temp));
      }else{
        likVal += AT[i]*log(pa) + AN[i]*log(1-pa) + BT[i]*log(pb) + BN[i]*log(1-pb) + temp;
      }
    }else{
      if ((pa==0 && AT[i]==0) || (pa==1 && AN[i]==0)){
        likVal += BT[i]*log(pb) + BN[i]*log(1-pb);
      }else if ((pa==0 && BT[i]==0) || (pa==1 && BN[i]==0)){
        likVal += AT[i]*log(pb) + AN[i]*log(1-pb);
      }else if ((pb==0 && AT[i]==0) || (pb==1 && AN[i]==0)){
        likVal += BT[i]*log(pa) + BN[i]*log(1-pa);
      }else if ((pb==0 && BT[i]==0) || (pb==1 && BN[i]==0)){
        likVal += AT[i]*log(pa) + AN[i]*log(1-pa);
      }else{
        likVal += log(0);
      }
    }
  }
  lik[0] = likVal;
  UNPROTECT(1);
  return(likS);
}
 
SEXP LikH(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP pS) {
  double *AT = REAL(ATS);
  double *BT = REAL(BTS);
  double *AN = REAL(ANS);
  double *BN = REAL(BNS);
  double *p = REAL(pS);
  double pa = p[0];
  double pb = p[1];
  long long len = length(ATS);
  SEXP likh;
  PROTECT(likh = allocVector(REALSXP, 2));
  double *likhPtr = REAL(likh);

  SEXP lik;
  PROTECT(lik = Lik(ATS, BTS, ANS, BNS, pS));
  likhPtr[0] = REAL(lik)[0];

  if (pa*(1-pa)*pb*(1-pb) ==0){
    likhPtr[1] = 0;
  }else{  
    // double H = 0;
    double ta = log(pa/(1-pa));
    double tb = log(pb/(1-pb));
    // double h1a = pa; // exp(ta)/(1+exp(ta));
    // double h1b = pb; // exp(tb)/(1+exp(tb));
    double h2a = pa*(1-pa); // exp(ta)/pow(1+exp(ta),2);
    double h2b = pb*(1-pb); // exp(tb)/pow(1+exp(tb),2);
    double l2a, l2ab, l2b, f12; // f12 = f1/f2
    l2a = l2ab = l2b = 0;
    for (long long i=0; i<len; i++){
      f12 = exp((AT[i]-BT[i])*(ta-tb) - (AT[i]+AN[i]-BT[i]-BN[i])*log((1-pb)/(1-pa)));
      // f1 = exp(AT[i]*ta - AN[i]*log(1+exp(ta)) + BT[i]*tb - BN[i]*log(1+exp(tb)));
      // f2 = exp(BT[i]*ta - BN[i]*log(1+exp(ta)) + AT[i]*tb - AN[i]*log(1+exp(tb)));
      l2a += pow((AT[i]-BT[i]-(AT[i]+AN[i]-BT[i]-BN[i])*pa),2)/(1+1/f12)/(f12+1) - h2a*((AT[i]+AN[i])/(1+1/f12)+(BT[i]+BN[i])/(1+f12));
      l2ab += (AT[i]-BT[i]-(AT[i]+AN[i]-BT[i]-BN[i])*pa)*(BT[i]-AT[i]-(BT[i]+BN[i]-AT[i]-AN[i])*pb)/(1+1/f12)/(1+f12);
      l2b += pow((BT[i]-AT[i]-(BT[i]+BN[i]-AT[i]-AN[i])*pb),2)/(1+1/f12)/(1+f12) - h2b*((BT[i]+BN[i])/(1+1/f12)+(AT[i]+AN[i])/(1+f12));
    }
    // printf("l2a: %f, l2b %f, l2ab %f\n", l2a, l2b, l2ab);
    likhPtr[1] = log(l2a*l2b - pow(l2ab,2));
  }
  UNPROTECT(2);
  return(likh);
}

SEXP SubSeq(SEXP ATS, long long begin, long long end) {
  // counting from 0, return ATS[begin, end)
  double *AT = REAL(ATS);
  SEXP newATS;
  PROTECT(newATS = allocVector(REALSXP, end-begin));
  double *newAT = REAL(newATS);
  for (long long i=0; i<(end-begin); i++){
    newAT[i] = AT[begin+i];
  }
  UNPROTECT(1);
  // double *test = REAL(newATS);
  
  return(newATS);
}

SEXP SubSeq2(SEXP ATS, long long begin, long long end) {
  // return the complement of ATS[begin, end)
  double *AT = REAL(ATS);
  long long len = length(ATS);
  SEXP newATS;
  PROTECT(newATS = allocVector(REALSXP, len-end+begin));
  double *newAT = REAL(newATS);
  for (long long i=0; i<begin; i++){
    newAT[i] = AT[i];
  }
  for (long long i=begin; i<len-end+begin; i++){
    newAT[i] = AT[end+i-begin];
  }
  UNPROTECT(1);
  // double *test = REAL(newATS);  
  return(newATS);
}

SEXP ScanStatNewCompBinom2dEMC(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP errorS, SEXP maxIterS, SEXP pS, SEXP gridCurS, SEXP maxWinS) {
  /* double *AT = REAL(ATS); */
  /* double *BT = REAL(BTS); */
  /* double *AN = REAL(ANS); */
  /* double *BN = REAL(BNS); */
  //  long long ATlen = length(ATS);
  SEXP p0;
  PROTECT(p0 = GetP(ATS, BTS, ANS, BNS, errorS, maxIterS, pS));
  long long maxWin = floor(REAL(maxWinS)[0]);
  double *gridCur = REAL(gridCurS);
  
  long long gridCurLen = length(gridCurS);
  long long gridCurMaxInd = gridCurLen - 1;
  long long rawi, rawj, i, j, jMax, bestWinI, bestWinJ;
  double l0, bestWinR, Rij;
  SEXP pij, pOut, newRes; //, ATij, BTij, ANij, BNij, ATOut, BTOut, ANOut, BNOut;
  SEXP l0S;
  PROTECT(l0S = Lik(ATS, BTS, ANS, BNS, p0));
  l0 = REAL(l0S)[0];
  UNPROTECT(2);  
  PROTECT(newRes = allocMatrix(REALSXP, gridCurMaxInd, 3));
  double *newResPtr = REAL(newRes);
  int newIter = 1;
  
  for (rawi=0; rawi<gridCurMaxInd; rawi++){
    jMax = rawi + maxWin;
    if (jMax > gridCurMaxInd)  jMax = gridCurMaxInd;
    bestWinI = gridCur[rawi];
    bestWinJ = gridCur[jMax];
    bestWinR = 0.0;
    newIter = 1;
    for (rawj=rawi+1; rawj<=jMax; rawj++){
      if (rawj-rawi == gridCurMaxInd) break;
      i = gridCur[rawi];
      j = gridCur[rawj]-1;
      SEXP ATij; PROTECT(ATij = SubSeq(ATS,i,j));
      SEXP BTij; PROTECT(BTij = SubSeq(BTS,i,j));
      SEXP ANij; PROTECT(ANij = SubSeq(ANS,i,j));
      SEXP BNij; PROTECT(BNij = SubSeq(BNS,i,j));
      PROTECT(pij = GetP(ATij, BTij, ANij, BNij, errorS, maxIterS, pS));
      SEXP lijS; 
      PROTECT(lijS = Lik(ATij, BTij, ANij, BNij, pij));
      Rij = REAL(lijS)[0];
      UNPROTECT(6); // ATij, BTij, ANij, BNij, pij, lij(Orig)      
      SEXP ATOut; PROTECT( ATOut = SubSeq2(ATS,i,j));
      SEXP BTOut; PROTECT( BTOut = SubSeq2(BTS,i,j));
      SEXP ANOut; PROTECT( ANOut = SubSeq2(ANS,i,j));
      SEXP BNOut; PROTECT( BNOut = SubSeq2(BNS,i,j));
      PROTECT(pOut = GetP(ATOut, BTOut, ANOut, BNOut, errorS, maxIterS, pS));
      SEXP lOutS;
      PROTECT(lOutS = Lik(ATOut, BTOut, ANOut, BNOut, pOut));
      Rij += REAL(lOutS)[0];
      UNPROTECT(6); // ATOut, BTOut, ANOut, BNOut, pOut, lOut(Orig)      
      if (Rij > bestWinR || newIter == 1) {
        bestWinI = i;
        bestWinJ = j+1;
        bestWinR = Rij;
        // printf("i: %lld, j: %lld, lij:%f, lOut:%f, Rij%f\n", i,j, lij, lOut, Rij);
      }
      newIter = 0;
    }
    bestWinR = bestWinR - l0;
    // if (bestWinR < 0) bestWinR = 0;
    newResPtr[rawi] = bestWinI;
    newResPtr[rawi + gridCurMaxInd] = bestWinJ;
    newResPtr[rawi + 2*gridCurMaxInd] = bestWinR;
  }
  UNPROTECT(1);
  return(newRes);
}
      
SEXP ScanStatRefineCompBinom2dEMC(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP errorS, SEXP maxIterS, SEXP pS, SEXP gridCurS, SEXP idLS, SEXP idRS) {
  /* double *AT = REAL(ATS); */
  /* double *BT = REAL(BTS); */
  /* double *AN = REAL(ANS); */
  /* double *BN = REAL(BNS); */
  SEXP p0;
  PROTECT(p0 = GetP(ATS, BTS, ANS, BNS, errorS, maxIterS, pS));
  long long gridCurLen = length(gridCurS);
  long long gridCurMaxInd = gridCurLen-1;
  double *gridCur = REAL(gridCurS);
  double *idL = REAL(idLS);
  double *idR = REAL(idRS);
  long long rawi, rawj, i, j, nRows, bestWinI, bestWinJ, rCt;
  double l0, lij, lOut, bestWinR, Rij;
  SEXP pij, pOut; //, ATij, BTij, ANij, BNij, ATOut, BTOut, ANOut, BNOut;
  int newIter = 1;
  SEXP l0S;
  bestWinI = 0;
  bestWinJ = 0;
  bestWinR = 0;
  PROTECT(l0S = Lik(ATS, BTS, ANS, BNS, p0));
  l0 = REAL(l0S)[0];
  UNPROTECT(2);
  nRows = length(idLS);
  if (idL[nRows-1] == gridCurMaxInd) nRows = nRows-1;
  // printf("nRows: %lld, length(idRS): %lld\n", nRows, (long long) length(idRS));
  
  SEXP newRes;
  PROTECT(newRes = allocMatrix(REALSXP, nRows, 3));
  double * newResPtr = REAL(newRes);
  rCt = 0;
  for (rawi=idL[0]; rawi<=idL[nRows-1]; rawi++){
    newIter = 1;
    for (rawj=idR[0]; rawj<=idR[length(idRS)-1]; rawj++){
      while(rawj < rawi+1){
        rawj++;
      }
      // printf("rawi:%lld, rawj:%lld\n", rawi, rawj);
      if (rawj-rawi == (length(gridCurS)-1)) break;
      i = gridCur[rawi];
      j = gridCur[rawj]-1;
      SEXP ATij; PROTECT(ATij = SubSeq(ATS,i,j));
      SEXP BTij; PROTECT(BTij = SubSeq(BTS,i,j));
      SEXP ANij; PROTECT(ANij = SubSeq(ANS,i,j));
      SEXP BNij; PROTECT(BNij = SubSeq(BNS,i,j));
      PROTECT(pij = GetP(ATij, BTij, ANij, BNij, errorS, maxIterS, pS));
      SEXP lijS; PROTECT(lijS = Lik(ATij, BTij, ANij, BNij, pij));
      lij = REAL(lijS)[0];
      UNPROTECT(6);
      Rij = lij;
      SEXP ATOut; PROTECT( ATOut = SubSeq2(ATS,i,j));
      SEXP BTOut; PROTECT( BTOut = SubSeq2(BTS,i,j));
      SEXP ANOut; PROTECT( ANOut = SubSeq2(ANS,i,j));
      SEXP BNOut; PROTECT( BNOut = SubSeq2(BNS,i,j));
      PROTECT(pOut = GetP(ATOut, BTOut, ANOut, BNOut, errorS, maxIterS, pS));
      SEXP lOutS; PROTECT(lOutS = Lik(ATOut, BTOut, ANOut, BNOut, pOut));
      lOut = REAL(lOutS)[0];
      UNPROTECT(6);
      Rij += lOut;
      // printf("Rij:%f\n", Rij);
      if (Rij > bestWinR || newIter == 1) {
        bestWinI = i;
        bestWinJ = j+1;
        bestWinR = Rij;
        // printf("i:%lld, j:%lld, Rij:%f\n", i, j, Rij);
     }
      newIter = 0;
    }
    bestWinR = bestWinR - l0;
    // if (bestWinR < 0) bestWinR = 0;
    newResPtr[rCt] = bestWinI;
    newResPtr[rCt + nRows] = bestWinJ;
    newResPtr[rCt + 2*nRows] = bestWinR;
    rCt++;
  }
  UNPROTECT(1);
  return(newRes);
}
