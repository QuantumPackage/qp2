#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "tree_utils.h"

void int_to_bin_digit(int64_t in, int count, int* out)
{
    /* assert: count <= sizeof(int)*CHAR_BIT */
    unsigned int mask = 1U << (count-1);
    int i;
    for (i = 0; i < count; i++) {
        out[i] = (in & mask) ? 1 : 0;
        in <<= 1;
    }
}

#include <stdio.h>
#include <stdint.h>
#include <math.h>

double logbinom(double n, double k) {
    return lgamma(n+1)-lgamma(n-k+1)-lgamma(k+1);
}
double binom(double n, double k) {
    return exp(logbinom(n,k));
}

void getncsfs1(int *inpnsomo, int *inpms, int *outncsfs){
    int nsomo = *inpnsomo;
    int ms = *inpms;
    int nparcoupl = (nsomo + ms)/2;
    *outncsfs = binom((double)nsomo, (double)nparcoupl);
}

void getncsfs(int NSOMO, int MS, int *outncsfs){
    int nparcoupl = (NSOMO + MS)/2; // n_alpha
    int nparcouplp1 = ((NSOMO + MS)/2)+1; // n_alpha + 1
    double tmpndets=0.0;
    if(NSOMO == 0){
        (*outncsfs) = 1;
        return;
    }
    tmpndets = binom((double)NSOMO, (double)nparcoupl);
    (*outncsfs) = round(tmpndets - binom((double)NSOMO, (double)nparcouplp1));
}

#include <stdint.h>

void getBFIndexList(int NSOMO, int *BF1, int *IdxListBF1){
    int Iidx;
    int Jidx;
    int BFcopy[NSOMO];

    int dictidx[2];
    dictidx[0] = -1;
    dictidx[1] =  1;

    for(int i = 0; i < NSOMO; i++)
        BFcopy[i] = BF1[i];

    for(int i = 0; i < NSOMO; i++){
        Iidx = i;
        if(BFcopy[i] == 0){
            int countN1=0;
            for(int j = i+1; j < NSOMO; j++){
                Jidx = j;
                countN1 = countN1 + dictidx[BFcopy[j]];
                if(countN1 > 0){
                    break;
                }
            }
            if(countN1 <= 0){
                BFcopy[Iidx] = -1;
                IdxListBF1[Iidx] = Iidx;
            }
            else{
                BFcopy[Iidx] = -1;
                BFcopy[Jidx] = -1;
                IdxListBF1[Jidx] = Iidx;
                IdxListBF1[Iidx] = Jidx;
            }
        }
    }

}

void getIslands(int NSOMO, int *BF1, int *BF2, int *nislands, int *phasefactor){

    // Get BF ids
    int *IdxListBF1 = malloc(NSOMO * sizeof(int));
    int *IdxListBF2 = malloc(NSOMO * sizeof(int));

    getBFIndexList(NSOMO, BF1, IdxListBF1);
    getBFIndexList(NSOMO, BF2, IdxListBF2);

    int sumids = 0;
    int maxcount=0;
    *nislands = 0;
    *phasefactor = 1;

    int BF1copy[NSOMO];
    for(int i = 0; i < NSOMO; i++)
        BF1copy[i] = IdxListBF1[i];
    int BF2copy[NSOMO];
    for(int i = 0; i < NSOMO; i++)
        BF2copy[i] = IdxListBF2[i];

    for(int i = 0; i < NSOMO; i++){
        int thisId = i;
        int nextId = BF1copy[i];
        maxcount = 0;
        while(BF1copy[thisId] != -1 && maxcount < 20){
            if(maxcount==0) *nislands += 1;
            if(maxcount==19) *nislands -= 1;

            maxcount++;

            // First the bra
            nextId = BF1copy[thisId];
            BF1copy[thisId] = -1;
            BF1copy[nextId] = -1;

            // Get the phase factor bra
            if(nextId < thisId) *phasefactor *= -1;

            // Then the ket
            thisId = BF2copy[nextId];
            BF2copy[thisId] = -1;
            BF2copy[nextId] = -1;

            // Get the phase factor bra
            if(nextId < thisId) *phasefactor *= -1;

        }
        
        for(int j=0;j<NSOMO;j++)
            sumids += BF1copy[j];
        if(sumids == -1*NSOMO) break;
        sumids = 0;
    }

    // Garbage collection
    free(IdxListBF1);
    free(IdxListBF2);

}

void getOverlapMatrix(int64_t Isomo, int64_t MS, double **overlapMatrixptr, int *rows, int *cols, int *NSOMOout){

    int NBF = 0;
    int NSOMO = 0;

    Tree bftree = (Tree){  .rootNode = NULL, .NBF = -1 };
    bftree.rootNode = malloc(sizeof(Node));
    (*bftree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    generateAllBFs(Isomo, MS, &bftree, &NBF, &NSOMO);

    *NSOMOout = NSOMO;

    // Initialize overlap matrix
    (*overlapMatrixptr) = malloc(NBF*NBF*sizeof(double));
    (*rows) = NBF;
    (*cols) = NBF;

    double *overlapMatrix = (*overlapMatrixptr);

    //// initialize Matrix
    //for(int i = 0; i < NBF; i++)
    //    for(int j = 0; j < NBF; j++)
    //        overlapMatrix[i*NBF + j] = 0.0;

    int addI = 0;
    int addJ = 0;
    int *BF1 = malloc(MAX_SOMO * sizeof(int));
    int *BF2 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF1 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF2 = malloc(MAX_SOMO * sizeof(int));

    int g = 0;
    g = (NSOMO - MS)/2;

    int nislands; // Note that nislands < g always
    int phasefactor;

    int dictPhase[2];

    dictPhase[0] = 1;
    dictPhase[1] =-1;


    // Set block elements
    for(int i = 0; i < NBF; i++){
        addI = i;
        getIthBFDriver(&bftree, NSOMO, addI, BF1);
        getBFIndexList(NSOMO, BF1, IdxListBF1);

        for(int j = 0; j < NBF; j++){
            addJ = j;
            getIthBFDriver(&bftree, NSOMO, addJ, BF2);
            getBFIndexList(NSOMO, BF2, IdxListBF2);

            // Get the i and r factors
            getIslands(NSOMO, BF1, BF2, &nislands, &phasefactor);

            overlapMatrix[i*NBF + j] = 1.0*phasefactor / (1 << (g - nislands));
        }
    }

    // Garbage collection
    free(BF1);
    free(IdxListBF1);
    free(BF2);
    free(IdxListBF2);

}


void getOverlapMatrix_withDet(double *bftodetmatrix, int rowsbftodetI, int colsbftodetI, int64_t Isomo, int64_t MS, double **overlapMatrixptr, int *rows, int *cols, int *NSOMOout){

    int NBF = 0;
    int NSOMO = 0;

    Tree bftree = (Tree){  .rootNode = NULL, .NBF = -1 };
    bftree.rootNode = malloc(sizeof(Node));
    (*bftree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    generateAllBFs(Isomo, MS, &bftree, &NBF, &NSOMO);

    (*NSOMOout) = NSOMO;

    // Initialize overlap matrix
    (*overlapMatrixptr) = malloc(NBF*NBF*sizeof(double));
    (*rows) = NBF;
    (*cols) = NBF;

    int transA=false;
    int transB=true;
    callBlasMatxMat(bftodetmatrix, rowsbftodetI, colsbftodetI, bftodetmatrix, rowsbftodetI, colsbftodetI, (*overlapMatrixptr), transA, transB);
}

void getSetBits(int64_t n, int *nsetbits){
    int count = 0;
    while(n){
        count += n & 1;
        n >>= 1;
    }
    *nsetbits = count;
}

void generateAllBFs(int64_t Isomo, int64_t MS, Tree *bftree, int *NBF, int *NSOMO){
    getSetBits(Isomo, NSOMO);
    buildTreeDriver(bftree, *NSOMO, MS, NBF);
}

//void ortho_qr_csf(double *overlapMatrix, int lda, double *orthoMatrix, int rows, int cols);


//void gramSchmidt_qp(double *overlapMatrix, int rows, int cols, double *orthoMatrix){
//  int i,j;
//  //for(j=0;j<cols;++j){
//  //  for(i=0;i<rows;++i){
//  //    printf(" %3.2f ",overlapMatrix[j*rows + i]);
//  //  }
//  //  printf("\n");
//  //}
//  // Call the function ortho_qr from qp
//  ortho_qr_csf(overlapMatrix, rows, orthoMatrix, rows, cols);
//  //for(j=0;j<cols;++j){
//  //  for(i=0;i<rows;++i){
//  //    printf(" %3.2f ",orthoMatrix[j*rows + i]);
//  //  }
//  //  printf("\n");
//  //}
//}

void gramSchmidt(double *overlapMatrix, int rows, int cols, double *orthoMatrix){

    // vector
    double norm = 0.0;
    double scalarprod = 0.0;
    orthoMatrix[(rows-1)*cols + cols-1] = 1.0;
    for(int i = cols-2; i > -1; i--){ orthoMatrix[(rows-1)*cols + i] = 0.0; }

    // Gram-Schmidt loop
    for(int i = rows-2; i > -1; i--){
        for(int k = cols-1; k > -1; k--){ orthoMatrix[(i)*cols + k] = 0.0; }
        orthoMatrix[i*cols + i] = 1.0;

        // orthogonalization
        for(int j = rows-1; j > i; j--){
            // calculate scalar product
            scalarprod = 0.0;
            for(int k = cols-1;k>=j;k--){
                scalarprod += orthoMatrix[j*cols + k] * overlapMatrix[i*cols + k];
            }
            for(int k = cols-1; k >= j; k--){
                orthoMatrix[i*cols + k] -= scalarprod * orthoMatrix[j*cols + k];
            }
        }

        // Normalization
        norm = 0.0;
        for(int j = rows-1; j >= i; j--){
            for(int k=cols-1; k >= i; k--)
                norm += orthoMatrix[i*cols + j]*orthoMatrix[i*cols + k]*overlapMatrix[j*cols+k];
        }
        norm = sqrt(norm);
        for(int j = rows-1; j >= i; j--){
            orthoMatrix[i*cols + j] /= norm;
        }

    }

}

void get_phase_cfg_to_qp_inpList(int *inpdet, int NSOMO, int *phaseout){
    int nbetas=0;
    (*phaseout) = 1;
    for(int i=0;i<NSOMO;i++){
        if(inpdet[i] == 0)
            (*phaseout) *= nbetas % 2 == 0 ? 1:-1;
        else
            nbetas += 1;
    }
    return;
}

void get_phase_cfg_to_qp_inpInt(int inpdet, double *phaseout){
    int nbetas=0;
    (*phaseout) = 1.0;
    int count=0;
    int mask=0;
    while(inpdet > 0){
        mask = (1<<count);
        if(__builtin_popcount(inpdet & mask)==1){
            (*phaseout) *= nbetas % 2 == 0 ? 1.0:-1.0;
            inpdet = inpdet ^ mask;
        }
        else nbetas += 1;
        count += 1;
    }
    //(*phaseout) = 1.0;
    return;
}

void convertCSFtoDetBasis(int64_t Isomo, int MS, int rowsmax, int colsmax, double *csftodetmatrix){

    double *overlapMatrixI;
    double *orthoMatrixI;
    double *bftodetmatrixI;
    double *csftodetmatrixI;
    int NSOMO=0;

    /***********************************
                 Get Overlap
    ************************************/
    // Fill matrix

    int rowsbftodetI, colsbftodetI;

    /***********************************
         Get BFtoDeterminant Matrix
    ************************************/


    //printf(" --- In convet ----\n");
    convertBFtoDetBasis(Isomo, MS, &bftodetmatrixI, &rowsbftodetI, &colsbftodetI);
    //printf(" --- done bf det basis ---- row=%d col=%d\n",rowsbftodetI,colsbftodetI);

    //printRealMatrix(bftodetmatrixI,rowsbftodetI,colsbftodetI);

    int rowsI = 0;
    int colsI = 0;

    //getOverlapMatrix(Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);
    getOverlapMatrix_withDet(bftodetmatrixI, rowsbftodetI, colsbftodetI, Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);

    //printf("Overlap matrix\n");
    //printRealMatrix(overlapMatrixI,rowsI,colsI);

    /***********************************
         Get Orthonormalization Matrix
    ************************************/

    orthoMatrixI = malloc(rowsI*colsI*sizeof(double));

    gramSchmidt(overlapMatrixI, rowsI, colsI, orthoMatrixI);

    //printf("Ortho matrix\n");
    //printRealMatrix(orthoMatrixI,rowsI,colsI);

    /***********************************
         Get Final CSF to Det Matrix
    ************************************/
    // First transform matrix using BLAS
    //double *bfIApqIJ = malloc(rowsbftodetI*colsbftodetI*sizeof(double));
    double *tmpcsftodet = malloc(rowsI*colsbftodetI*sizeof(double));

    int transA=false;
    int transB=false;
    double phaseAll = -1.0;
    callBlasMatxMat(orthoMatrixI, rowsI, colsI, bftodetmatrixI, rowsbftodetI, colsbftodetI, tmpcsftodet, transA, transB);
    for(int i=0;i<rowsI;i++){
        phaseAll = 1.0;
        //for(int j=0;j<colsbftodetI;j++){
        //    if(tmpcsftodet[i*colsbftodetI + j] > 0.0) phaseAll = 1.0;
        //}
        for(int j=0;j<colsbftodetI;j++){
            csftodetmatrix[j*rowsI + i] = tmpcsftodet[i*colsbftodetI + j]*phaseAll;
        }
    }

    // Garbage collection
    if(rowsI + colsI > 0) free(overlapMatrixI);
    if(rowsI + colsI > 0) free(orthoMatrixI);
    if(rowsbftodetI + colsbftodetI > 0) free(bftodetmatrixI);
    if(rowsI + colsbftodetI > 0) free(tmpcsftodet);
}

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0')

int applyRemoveShftAddSOMOVMO(int idet, int p, int q, int *phase){
    // CSF: 1 1 1 1 0 1
    // DET: 1 0 1 0   1
    //        |     |
    //        p     q
    //        p = 4
    //        q = 1
    //
    //          result
    //
    // CSF: 1 0 1 1 1 1
    // DET: 1   1 0 0 1
    // maskp:
    //      0   1 1 1 1
    // maskq:
    //      0   0 0 0 1
    // maskpxq:
    //      0   1 1 1 0
    // maskqxqi:
    //      1   0 0 0 1
    int maskp  = (1UL << p)-1;
    int maskq  = (1UL << q)-1;
    int maskpxq = (maskp ^ maskq);
    int maskpxqi = ~(maskp ^ maskq);

    // Step 1: remove
    // clear bits from p
    int outdet = idet;
    int occatp = __builtin_popcount(idet & (1UL << (p-1)));
    // remove the bit at p
    outdet &= ~(1UL << (p-1));

    // Step 2: shift
    if(q > p){
        // start with q

        // calculate the phase
        int na, nb;
        int tmpdet = outdet & (maskpxq);
        na = __builtin_popcount(tmpdet);
        nb = __builtin_popcount(maskpxq) - na;
        //int nfermions = occatp == 0 ? nb : na;
        int nfermions = na+nb;
        (*phase) = nfermions % 2 == 0 ? 1 : -1;

        int tmpdetq1 = outdet & maskpxq;
        int tmpdetq2 = outdet & maskpxqi;
        tmpdetq1 = tmpdetq1 >> 1;
        outdet = tmpdetq1 | tmpdetq2;
        // put electron at q
        outdet = occatp == 0 ? outdet : outdet | (1UL<<(q-1));
    }
    else{
        // shift bit to right
        maskpxq = maskpxq >> 1;
        maskpxqi = ~(maskpxq);

        // calculate the phase
        int na, nb;
        int tmpdet = outdet & (maskpxq);
        na = __builtin_popcount(tmpdet);
        nb = __builtin_popcount(maskpxq) - na;
        //int nfermions = occatp == 0 ? nb : na;
        int nfermions = na+nb;
        (*phase) = nfermions % 2 == 0 ? 1 : -1;

        // start with p
        // shift middle electrons to right
        int tmpdetp1 = outdet & maskpxq;
        int tmpdetp2 = outdet & maskpxqi;
        tmpdetp1 = tmpdetp1 << 1;
        outdet = tmpdetp1 | tmpdetp2;
        // put electron at q
        outdet = occatp == 0 ? outdet : outdet | (1UL<<(q-1));
    }

    // Done
    return(outdet);
}

int applyRemoveShftAddDOMOSOMO(int idet, int p, int q, int *phase){
    // CSF: 1 2 1 1 1 1 1 1 1 1
    // DET: 1   0 0 1 1 0 0 1 0
    //          |       |
    //          p       q
    //
    //          result
    //
    // CSF: 1 1 1 1 1 1 2 1 1 1
    // DET: 1 0 0 0 1 1   0 1 0
    // maskp:
    //      0   1 1 1 1 1 1 1 1
    // maskq:
    //      0 0 0 0 0 0 1 1 1 1
    int maskp  = (1UL << p)-1;
    int maskq  = (1UL << q)-1;
    int maskpxq = (maskp ^ maskq);
    int maskpxqi = ~(maskp ^ maskq);

    // Step 1: remove
    // clear bits from q
    int outdet = idet;
    int occatq = __builtin_popcount(idet & (1UL << (q-1)));
    outdet &= ~(1UL << (q-1));

    // Step 2: shift
    if(q > p){
        // start with q

        // shift mask between p and q
        maskpxq = maskpxq >> 1;
        maskpxqi = ~(maskpxq);
        // calculate the phase
        int na, nb;
        int tmpdet = outdet & (maskpxq);
        na = __builtin_popcount(tmpdet);
        nb = __builtin_popcount(maskpxq) - na;
        // spin obb to that at q is moving
        //int nfermions = occatq == 0 ? na : nb;
        int nfermions = na + nb + 1;
        (*phase) = nfermions % 2 == 0 ? 1 : -1;

        int tmpdetq1 = outdet & maskpxq;
        int tmpdetq2 = outdet & maskpxqi;
        tmpdetq1 = tmpdetq1 << 1;
        outdet = tmpdetq1 | tmpdetq2;

        // Step 3: Add bit at p + 1
        outdet = occatq == 1 ? outdet | (1UL<<(p-1)) : outdet;
    }
    else{

        // calculate the phase
        int na, nb;
        int tmpdet = outdet & (maskpxq);
        na = __builtin_popcount(tmpdet);
        nb = __builtin_popcount(maskpxq) - na;
        // spin obb to that at q is moving
        //int nfermions = occatq == 0 ? na : nb;
        int nfermions = na + nb + 1;
        (*phase) = nfermions % 2 == 0 ? 1 : -1;

        // start with p
        // shift middle electrons to right
        int tmpdetp1 = outdet & maskpxq;
        int tmpdetp2 = outdet & maskpxqi;
        tmpdetp1 = tmpdetp1 >> 1;
        outdet = tmpdetp1 | tmpdetp2;

        // Step 3: Add bit at p
        outdet = occatq == 1 ? outdet | (1UL<<(p-1)) : outdet;
    }

    // Done
    return(outdet);
}

int applyRemoveShftSOMOSOMO(int idet, int p, int q, int *phase){
    // CSF: 1 1 1 1 1 1 1 1 1 1
    // DET: 1 1 0 0 1 1 0 0 1 0
    //        |         |
    //        p         q
    //
    //          result
    //
    // CSF: 1   1 1 1 1   1 1 1
    // DET: 1   0 0 1 1   0 1 0
    // maskp:
    //      0 1 1 1 1 1 1 1 1 1
    // maskq:
    //      0 0 0 0 0 0 0 1 1 1
    int maskp  = (1UL << p)-1;
    int maskq  = (1UL << q)-1;
    int maskpi =~maskp;
    int maskqi =~maskq;

    // Step 1: remove
    // clear bits from p and q
    int outdet = idet;
    outdet &= ~(1UL << (p-1));
    outdet &= ~(1UL << (q-1));

    // calculate the phase
    int occatp = idet & (1UL << (p-1));
    int na, nb;
    int tmpdet = outdet & (maskp ^ maskq);
    na = __builtin_popcount(tmpdet);
    nb = abs(p-q)-1 - na;
    //int nfermions = occatp == 0 ? nb : na;

    // Step 2: shift
    if(q > p){
        int nfermions = occatp == 0 ? na+nb : na+nb+1;
        (*phase) = nfermions % 2 == 0 ? 1 : -1;
        // start with q
        // shift everything left of q
        int tmpdetq1 = outdet & maskq;
        int tmpdetq2 = outdet & maskqi;
        tmpdetq2 = tmpdetq2 >> 1;
        outdet = tmpdetq1 | tmpdetq2;

        // shift everything left of p
        int tmpdetp1 = outdet & maskp;
        int tmpdetp2 = outdet & maskpi;
        tmpdetp2 = tmpdetp2 >> 1;
        outdet = tmpdetp1 | tmpdetp2;
    }
    else{
        int nfermions = occatp == 0 ? na+nb+1 : na+nb;
        (*phase) = nfermions % 2 == 0 ? 1 : -1;
        // start with p
        // shift everything left of p
        int tmpdetp1 = outdet & maskp;
        int tmpdetp2 = outdet & maskpi;
        tmpdetp2 = tmpdetp2 >> 1;
        outdet = tmpdetp1 | tmpdetp2;

        // shift everything left of q
        int tmpdetq1 = outdet & maskq;
        int tmpdetq2 = outdet & maskqi;
        tmpdetq2 = tmpdetq2 >> 1;
        outdet = tmpdetq1 | tmpdetq2;
    }

    // Done
    return(outdet);
}

unsigned int shftbit(int num, int p){
    unsigned int maskleft = ~(0 | ((1<<p)-1));
    unsigned int maskright = ((1<<(p-1))-1);
    int numleft = num & maskleft;
    int numright = num & maskright;
    numleft = numleft >> 1;
    return(numleft | numright);
};

int getphase(int num, int p, int q, int nmo){
    // CSF: 1 1 1 1 1 1 1 1 1 1
    // DET: 1 1 0 0 1 1 0 0 1 0
    //        |         |
    //        p         q
    //        |         |
    // CSF: 1 1 1 1 1 1 1 1 1 1
    // DET: 1 0 0 0 1 1 1 0 1 0
    //
    // maskleft:
    //      1 1 1 1 1 1 1 0 0 0
    // maskright:
    //      0 1 1 1 1 1 1 1 1 1
    int omax = p > q ? p : q;
    int omin = p > q ? q : p;
    unsigned int maskleft = ~(0 | ((1<<(omin-1))-1));
    unsigned int maskright = ((1<<(omax))-1);
    unsigned int maskmo = ((1<<nmo)-1);
    int numleft = num & maskleft;
    int numleftright = numleft & maskright;
    int nalpha = __builtin_popcount(numleftright & maskmo);
    int nbeta = omax-omin+1 - nalpha;
    int maskatp = (1<<(p-1));
    int nelecalphaatp = __builtin_popcount(num & maskatp);
    int maskatq = (1<<(q-1));
    int nelecalphaatq = __builtin_popcount(num & maskatq);
    int nfermions = nelecalphaatp == 0 ? nbeta : nalpha;
    int phase = (nfermions-1) % 2 == 0 ? 1 : -1;
    if(nelecalphaatp == nelecalphaatq) phase = 0.0;
    return(phase);
};


int getDOMOSOMOshift(int idet, int p, int q, int *phase){
    /*
      Idea:
      DOMO->SOMO example

      1 2 1 1 1
        p     q
      1 1 1 1 2

      p = 3
      q = 1

      in determinant representation: (0->beta,1->alpha)
      |I>   = 0 0 1 1
               |____|
               p    q

      |ret> = 0 1 0 1
      A shift of bit at q to pos after p.

    */

    int maskq = ~((1UL<<q)-1);
    int maskp = (1UL<<p)-1;
    int maskpq = ~(maskp & maskq);
    int bits_to_shft = (idet & maskq) & maskp;
    // shift bits by 1 index
    int shifted_bits = bits_to_shft >> 1;
    // Now combine with original det
    int detout = (idet & maskpq);
    // Zero out bits at q
    detout &= ~(1UL << (q-1));
    // Set the bit at p
    detout |=  (1UL << (p-1));
    // Add the shifted bits
    detout |= shifted_bits;

    // Now calcaulate the phase
    // Find the type of bit at q
    int occatq = idet & (1UL << (q-1));
    // calculate number of alpha and beta spins
    int na = __builtin_popcount(shifted_bits);
    int nb = p - q - na;
    // Find the number of fermions to pass
    int nfermions = occatq == 0 ? na : nb;
    (*phase) = nfermions % 2 == 0 ? 1 : -1;
    return(detout);
}

void calcMEdetpair(int *detlistI, int *detlistJ, int orbI, int orbJ, int Isomo, int Jsomo, int ndetI, int ndetJ, int NMO, double *matelemdetbasis){

    // Calculation of phase
    // The following convention is used
    // <J|a^{\dagger}_q a_p | I>
    //
    // The phase is calculated
    // assuming all alpha electrons
    // are on the left and all beta
    // electrons are on the RHS
    // of the alphas.


    int maskI;
    int nelecatI;
    unsigned int maskleft;
    unsigned int maskright;
    unsigned int psomo;
    unsigned int qsomo;


    // E(q,p) |I> = cqp |J>


    int p,q; // The two orbitals p is always > q.
    p = orbI >= orbJ ? orbI : orbJ;
    q = orbI >= orbJ ? orbJ : orbI;

    // Find the corresponding case
    // 1. NdetI > NdetJ  (SOMO -> SOMO)
    // 2. NdetI < NdetJ  (DOMO -> VMO)
    // 3. NdetI == NdetJ (SOMO -> VMO and DOMO -> SOMO)

    // Converting the above four cases into int:
    int case_type = abs(ndetI - ndetJ) == 0 ? 3 : (ndetI > ndetJ ? 1 : 2);

    switch (case_type){
        case 1:
            // SOMO -> SOMO
            // Find the orbital ids in model space
            maskleft  =  (0 | ((1<<(p))-1));
            maskright =  (0 | ((1<<(q))-1));
            psomo = __builtin_popcount(Isomo & maskleft);
            qsomo = q == 1 ? 1 : __builtin_popcount(Isomo & maskright);
            p = psomo >= qsomo ? psomo : qsomo;
            q = psomo >= qsomo ? qsomo : psomo;

            for(int i=0;i<ndetI;i++){
                int idet = detlistI[i];
                int phase = getphase(idet,orbI,orbJ,NMO);
                // Shift bits for
                idet = shftbit(shftbit(detlistI[i],q),p-1);
                for(int j=0;j<ndetJ;j++){
                    int jdet = (detlistJ[j]);
                    if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                }
            }
            break;
        case 2:
            // DOMO -> VMO
            // Find the orbital ids in model space
            maskleft = (0 | ((1<<(p))-1));
            maskright =(0 | ((1<<(q))-1));
            psomo = __builtin_popcount(Jsomo & maskleft);
            qsomo = q == 1 ? 1 : __builtin_popcount(Jsomo & maskright);
            p = psomo >= qsomo ? psomo : qsomo;
            q = psomo >= qsomo ? qsomo : psomo;

            for(int i=0;i<ndetI;i++){
                // Get phase
                int idet = detlistI[i];
                for(int j=0;j<ndetJ;j++){
                    int jdet = (detlistJ[j]);
                    // Calculate phase
                    int phase = 1*getphase(jdet,p,q,NMO);
                    // Shift bits for I
                    jdet = shftbit(shftbit(detlistJ[j],q),p-1);
                    if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                }
            }
            break;
        case 3:
            // (SOMO -> VMO or DOMO -> SOMO)
            // if Isomo[p] == 1 => SOMO -> VMO
            // if Isomo[p] == 0 => DOMO -> SOMO
            // Find the orbital ids in model space
            maskleft = ((1<<(p))-1);
            maskright =((1<<(q))-1);
            psomo = __builtin_popcount(Isomo & maskleft);
            //qsomo = q == 1 ? 1 : __builtin_popcount(Isomo & maskright);
            qsomo = __builtin_popcount(Isomo & maskright);
            p = psomo >= qsomo ? psomo : qsomo;
            q = psomo >= qsomo ? qsomo : psomo;


            int noccorbI = (Isomo & (1<<(orbI-1)));
            switch (noccorbI){
                case 0:
                    // Case: DOMO -> SOMO
                    break;
                case 1:
                    // Case: SOMO -> VMO
                    break;
                default:
                    printf("Something is wrong in calcMEdetpair\n");
                    break;
            }

            int tmpidet;

            for(int i=0;i<ndetI;i++){
                // Get phase
                int idet = detlistI[i];
                int nelecalphaatp = (Isomo & (1<<(orbI-1)));
                // Idea:
                // if DOMO -> SOMO
                //
                // I =
                //   2  1 1 1 1
                // (10) 0 0 1 1
                //
                //     |
                //    \ /
                //     .
                //  0 0 0 1 1
                //
                // J =
                // 1 1 1 1  2
                // 0 0 1 1 (10)
                //
                if(nelecalphaatp == 0){
                    // Case: DOMO -> SOMO
                    tmpidet = idet;
                    int nelecalphaatq = (idet & (1<<(orbJ-1)));
                    if(nelecalphaatq==0) tmpidet = tmpidet ^ (1<<(orbI-1));
                    else                 tmpidet = tmpidet ^ (0);
                    idet = shftbit(idet,q);
                }
                else{
                    tmpidet = idet;
                    idet = shftbit(idet,p);
                }

                // Calculate phase
                int phase = 1*getphase(tmpidet,orbI,orbJ,NMO);
                for(int j=0;j<ndetJ;j++){
                    int jdet;
                    if(nelecalphaatp == 0) jdet = shftbit(detlistJ[j],p);
                    else                   jdet = shftbit(detlistJ[j],q);
                    if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                }
            }

            break;
        default:
            printf("Something is wrong in calc ME\n");
            break;
    } // end select

}

void calcMEdetpairGeneral(int *detlistI, int *detlistJ, int orbI, int orbJ, int Isomo, int Jsomo, int ndetI, int ndetJ, int NMO, double *matelemdetbasis){

    // Calculation of phase
    // The following convention is used
    // <J|a^{\dagger}_q a_p | I>
    //
    // The phase is calculated
    // assuming all alpha electrons
    // are on the left and all beta
    // electrons are on the RHS
    // of the alphas.

    // There are three possibilities
    // which need to be separated
    // CASE 1. p > q
    // CASE 2. p < q
    // CASE 3. p == q

    int maskI;
    int nelecatI;
    int noccorbI;
    double phaseI=1.0;
    double phaseJ=1.0;
    unsigned int maskleft;
    unsigned int maskright;
    unsigned int psomo;
    unsigned int qsomo;

    int p,q; // The two orbitals p is always > q.

    if(orbI > orbJ){
        // CASE 1 : orbI > orbJ
        p = orbI;
        q = orbJ;

        // Find the corresponding sub case
        // 1. NdetI > NdetJ  (SOMO -> SOMO)
        // 2. NdetI < NdetJ  (DOMO -> VMO)
        // 3. NdetI == NdetJ (SOMO -> VMO and DOMO -> SOMO)

        // Converting the above four cases into int:
        int case_type = abs(ndetI - ndetJ) == 0 ? 3 : (ndetI > ndetJ ? 1 : 2);
        p = orbI;
        q = orbJ;

        switch (case_type){
            case 1:
                // SOMO -> SOMO
                // Find the orbital ids in model space
                maskleft  =  (0 | ((1<<(p))-1));
                maskright =  (0 | ((1<<(q))-1));
                psomo = __builtin_popcount(Isomo & maskleft);
                qsomo = __builtin_popcount(Isomo & maskright); // q has to be atleast 1
                p = psomo;
                q = qsomo;

                for(int i=0;i<ndetI;i++){
                    int idet = detlistI[i];
                    int phase=1;
                    // Apply remove and shft on Isomo
                    idet = applyRemoveShftSOMOSOMO(idet, p, q, &phase);
                    //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                    for(int j=0;j<ndetJ;j++){
                        int jdet = (detlistJ[j]);
                        //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                        if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                    }
                }
                break;
            case 2:
                // DOMO -> VMO
                // Find the orbital ids in model space
                // As seen in Jsomo
                // Here we apply a^{\dagger}_p a_q |J>
                maskleft = (0 | ((1<<(p))-1));
                maskright =(0 | ((1<<(q))-1));
                psomo = __builtin_popcount(Jsomo & maskleft);
                qsomo = __builtin_popcount(Jsomo & maskright); // q has to be atleast 1
                p = psomo;
                q = qsomo;

                for(int i=0;i<ndetI;i++){
                    // Get phase
                    int idet = detlistI[i];
                    //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                    for(int j=0;j<ndetJ;j++){
                        int jdet = (detlistJ[j]);
                        //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                        // Calculate phase
                        int phase=1;
                        // Apply remove and shift on Jdet (orbital ids are inverted)
                        jdet = applyRemoveShftSOMOSOMO(jdet, q, p, &phase);
                        if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                    }
                }
                break;
            case 3:
                // (SOMO -> VMO or DOMO -> SOMO)
                noccorbI = __builtin_popcount(Isomo & (1<<(orbI-1)));

                switch (noccorbI){
                    case 0:
                        // Case: DOMO -> SOMO
                        // Find the orbital ids in model space
                        // Ex:
                        //      2 1 1 1 1
                        //      p     q
                        //      1 1 1 2 1
                        // p = 4
                        // q = 2
                        // p is from Jsomo
                        // q is from Isomo
                        maskleft = ((1<<(p))-1);
                        maskright =((1<<(q))-1);
                        psomo = __builtin_popcount(Jsomo & maskleft);
                        qsomo = __builtin_popcount(Isomo & maskright);
                        p = psomo;
                        q = qsomo;

                        for(int i=0;i<ndetI;i++){
                            int idet = detlistI[i];
                            //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                            int phase=1;
                            // Apply remove and shft on Isomo
                            idet = applyRemoveShftAddDOMOSOMO(idet, p, q, &phase);
                            for(int j=0;j<ndetJ;j++){
                                int jdet = (detlistJ[j]);
                                //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                                if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                            }
                        }
                        break;
                    case 1:
                        // Case: SOMO -> VMO
                        // Find the orbital ids in model space
                        // Ex:
                        //      1 1 1 0 1
                        //      p     q
                        //      0 1 1 1 1
                        // p = 4
                        // q = 1
                        // p is from Isomo
                        // q is from Jsomo
                        maskleft = ((1<<(p))-1);
                        maskright =((1<<(q))-1);
                        psomo = __builtin_popcount(Isomo & maskleft);
                        qsomo = __builtin_popcount(Jsomo & maskright);
                        p = psomo;
                        q = qsomo;

                        for(int i=0;i<ndetI;i++){
                            int idet = detlistI[i];
                            //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                            int phase=1;
                            // Apply remove and shft on Isomo
                            idet = applyRemoveShftAddSOMOVMO(idet, p, q, &phase);
                            for(int j=0;j<ndetJ;j++){
                                int jdet = (detlistJ[j]);
                                //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                                if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                            }
                        }
                        break;
                    default:
                        printf("Something is wrong in calcMEdetpair\n");
                        break;
                }
                break;
            default:
                printf("Something is wrong in calc ME\n");
                break;
        } // end select

    } // end orbI > orbJ
    else if(orbI < orbJ){
        // CASE 2 orbI < orbJ
        p = orbI;
        q = orbJ;
        // Find the corresponding sub case
        // 1. NdetI > NdetJ  (SOMO -> SOMO)
        // 2. NdetI < NdetJ  (DOMO -> VMO)
        // 3. NdetI == NdetJ (SOMO -> VMO and DOMO -> SOMO)

        // Converting the above four cases into int:
        int case_type = abs(ndetI - ndetJ) == 0 ? 3 : (ndetI > ndetJ ? 1 : 2);

        switch (case_type){
            case 1:
                // SOMO -> SOMO
                // Find the orbital ids in model space
                maskleft  =  (0 | ((1<<(p))-1));
                maskright =  (0 | ((1<<(q))-1));
                psomo = __builtin_popcount(Isomo & maskleft);
                qsomo = __builtin_popcount(Isomo & maskright); // q has to be atleast 1
                p = psomo;
                q = qsomo;

                for(int i=0;i<ndetI;i++){
                    int idet = detlistI[i];
                    //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                    int phase=1;
                    // Apply remove and shft on Isomo
                    idet = applyRemoveShftSOMOSOMO(idet, p, q, &phase);
                    for(int j=0;j<ndetJ;j++){
                        int jdet = (detlistJ[j]);
                        //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                        if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                    }
                }
                break;
            case 2:
                // DOMO -> VMO
                // Find the orbital ids in model space
                // As seen in Jsomo
                // Here we apply a^{\dagger}_p a_q |J>
                maskleft = (0 | ((1<<(p))-1));
                maskright =(0 | ((1<<(q))-1));
                psomo = __builtin_popcount(Jsomo & maskleft);
                qsomo = __builtin_popcount(Jsomo & maskright); // q has to be atleast 1
                p = psomo;
                q = qsomo;

                for(int i=0;i<ndetI;i++){
                    // Get phase
                    int idet = detlistI[i];
                    //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                    for(int j=0;j<ndetJ;j++){
                        int jdet = (detlistJ[j]);
                        //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                        // Calculate phase
                        int phase=1;
                        // Apply remove and shift on Jdet (orbital ids are inverted)
                        jdet = applyRemoveShftSOMOSOMO(jdet, q, p, &phase);
                        if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                    }
                }
                break;
            case 3:
                // (SOMO -> VMO or DOMO -> SOMO)
                // if Isomo[p] == 1 => SOMO -> VMO
                // if Isomo[p] == 0 => DOMO -> SOMO
                noccorbI = __builtin_popcount(Isomo & (1<<(orbI-1)));

                switch (noccorbI){
                    case 0:
                        // Case: DOMO -> SOMO
                        // Find the orbital ids in model space
                        // Ex:
                        //      1 1 1 2 1
                        //      q     p
                        //      2 1 1 1 1
                        // p = 1
                        // q = 4
                        // p is from Jsomo
                        // q is from Isomo
                        maskleft = ((1<<(p))-1);
                        maskright =((1<<(q))-1);
                        psomo = __builtin_popcount(Jsomo & maskleft);
                        qsomo = __builtin_popcount(Isomo & maskright);
                        p = psomo;
                        q = qsomo;

                        for(int i=0;i<ndetI;i++){
                            int idet = detlistI[i];
                            //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                            int phase=1;
                            // Apply remove and shft on Isomo
                            idet = applyRemoveShftAddDOMOSOMO(idet, p, q, &phase);
                            for(int j=0;j<ndetJ;j++){
                                int jdet = (detlistJ[j]);
                                //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                                if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                            }
                        }
                        break;
                    case 1:
                        // Case: SOMO -> VMO
                        // Find the orbital ids in model space
                        // Ex:
                        //      0 1 1 1 1
                        //      q     p
                        //      1 1 1 0 1
                        // p = 2
                        // q = 4
                        // p is from Isomo
                        // q is from Jsomo
                        maskleft = ((1<<(p))-1);
                        maskright =((1<<(q))-1);
                        psomo = __builtin_popcount(Isomo & maskleft);
                        qsomo = __builtin_popcount(Jsomo & maskright);
                        p = psomo;
                        q = qsomo;

                        for(int i=0;i<ndetI;i++){
                            int idet = detlistI[i];
                            //get_phase_cfg_to_qp_inpInt(detlistI[i], &phaseI);
                            int phase=1;
                            // Apply remove and shft on Isomo
                            idet = applyRemoveShftAddSOMOVMO(idet, p, q, &phase);
                            for(int j=0;j<ndetJ;j++){
                                int jdet = (detlistJ[j]);
                                //get_phase_cfg_to_qp_inpInt(detlistJ[j], &phaseJ);
                                if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                            }
                        }
                        break;
                    default:
                        printf("Something is wrong in calcMEdetpair\n");
                        break;
                }
                break;
            default:
                printf("Something is wrong in calc ME\n");
                break;
        } // end select
    } // end orbI  < orbJ
    else{
        // CASE 3 : orbI == orbJ

        // Three possibilities
        // orbI = VMO
        // orbI = SOMO
        // orbI = DOMO
        int noccorbI = __builtin_popcount((Isomo & (1<<(orbI-1))));
        switch (noccorbI){
            case 0:
                // Matrix is 0
                for(int i=0;i<ndetI;i++){
                    int idet = detlistI[i];
                    for(int j=0;j<ndetJ;j++){
                        int jdet = (detlistJ[j]);
                        if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 0.0;
                    }
                }
                break;
            case 1:
                // Matrix is Identity
                for(int i=0;i<ndetI;i++){
                    int idet = detlistI[i];
                    for(int j=0;j<ndetJ;j++){
                        int jdet = (detlistJ[j]);
                        if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0;
                    }
                }
                break;
            default:
                break;
        }
    } // end orbI == orbJ

    return;
}

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0')

void callcalcMEij(int Isomo, int Jsomo, int orbI, int orbJ, int MS, int NMO, double **ApqIJptr, int *rowsA, int *colsA){
    // Get dets for I
    int ndetI;
    int ndetJ;

    // Get detlist
    int NSOMOI=0;
    int NSOMOJ=0;
    getSetBits(Isomo, &NSOMOI);
    getSetBits(Jsomo, &NSOMOJ);

    Tree dettreeI = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettreeI.rootNode = malloc(sizeof(Node));
    (*dettreeI.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettreeI, Isomo, MS, &ndetI);


    Tree dettreeJ = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettreeJ.rootNode = malloc(sizeof(Node));
    (*dettreeJ.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettreeJ, Jsomo, MS, &ndetJ);

    int detlistI[ndetI];
    int detlistJ[ndetJ];
    for(int i=0;i<ndetI;i++)
        detlistI[i] = 0;
    for(int i=0;i<ndetJ;i++)
        detlistJ[i] = 0;

    // Get detlist
    getDetlistDriver(&dettreeI, NSOMOI, detlistI);
    getDetlistDriver(&dettreeJ, NSOMOJ, detlistJ);

    (*ApqIJptr) = malloc(ndetI*ndetJ*sizeof(double));
    (*rowsA) = ndetI;
    (*colsA) = ndetJ;

    double *matelemdetbasis = (*ApqIJptr);

    for(int i=0;i<ndetI;i++)
        for(int j=0;j<ndetJ;j++)
            matelemdetbasis[i*ndetJ + j]=0.0;

    // Calculate matrix elements in det basis
    //calcMEdetpair(detlistI, detlistJ, orbI, orbJ, Isomo, Jsomo, ndetI, ndetJ, NMO, matelemdetbasis);
    calcMEdetpairGeneral(detlistI, detlistJ, orbI, orbJ, Isomo, Jsomo, ndetI, ndetJ, NMO, matelemdetbasis);

}

void getbftodetfunction(Tree *dettree, int NSOMO, int MS, int *BF1, double *rowvec){
    int npairs = 1 << ((NSOMO - MS)/2);
    int idxp = 0;
    int idxq = 0;
    int *detslist = malloc(npairs*NSOMO*sizeof(int));
    double *phaselist = malloc(npairs*sizeof(double));
    for(int i=0;i<npairs;i++)
        phaselist[i] = 1.0;
    int shft = npairs;
    int donepq[NSOMO];
    double fac = 1.0;
    for(int i = 0; i < NSOMO; i++)
        donepq[i] = 0.0;
    for(int i=0;i<npairs;++i){
      for(int j=0;j<NSOMO;++j)
        detslist[i*NSOMO + j]=0;
    }

    for(int i = 0; i < NSOMO; i++){
        idxp = BF1[i];
        idxq = BF1[idxp];
        // Do one pair only once
        if(donepq[idxp] > 0.0 || donepq[idxq] > 0.0 || idxp == idxq) continue;
        fac *= 2.0;
        donepq[idxp] = 1.0;
        donepq[idxq] = 1.0;
        for(int j = 0; j < npairs; j = j + shft){
            for(int k = 0; k < shft/2; k++){
                detslist[(k+j)*NSOMO + idxp] = 1;
                detslist[(k+j)*NSOMO + idxq] = 0;
            }
            for(int k = shft/2; k < shft; k++){
                detslist[(k+j)*NSOMO + idxp] = 0;
                detslist[(k+j)*NSOMO + idxq] = 1;
                phaselist[k+j] *=-1;
            }
        }
        shft /= 2;
    }
    
    // Now get the addresses
    int inpdet[NSOMO];
    int phase_cfg_to_qp=1;
    int addr = -1;
    for(int i = 0; i < npairs; i++){
        for(int j = 0; j < NSOMO; j++) {
            inpdet[j] = detslist[i*NSOMO + j];
            //printf(" %d ",inpdet[j]);
        }
        //printf("\n");
        findAddofDetDriver(dettree, NSOMO, inpdet, &addr);
        //printf("(%d) - addr  = %d\n",i,addr);
        // Calculate the phase for cfg to QP2 conversion
        //get_phase_cfg_to_qp_inpList(inpdet, NSOMO, &phase_cfg_to_qp);
        //rowvec[addr] = 1.0 * phaselist[i]*phase_cfg_to_qp/sqrt(fac);
        rowvec[addr] = 1.0 * phaselist[i]/sqrt(fac);
        // Upon transformation from
        // SOMO to DET basis,
        // all dets have the same phase
        // Is this true ?
        //rowvec[addr] = 1.0/sqrt(fac);
    }

    free(detslist);
    free(phaselist);
}

void convertBFtoDetBasis(int64_t Isomo, int MS, double **bftodetmatrixptr, int *rows, int *cols){

    int NSOMO=0;
    //printf("before getSetBits Isomo=%ld, NSOMO=%ld\n",Isomo,NSOMO);
    getSetBits(Isomo, &NSOMO);
    //printf("Isomo=%ld, NSOMO=%ld\n",Isomo,NSOMO);
    int ndets = 0;
    int NBF = 0;
    //double dNSOMO = NSOMO*1.0;
    // MS = alpha_num - beta_num
    int nalpha = (NSOMO + MS)/2;
    //printf(" in convertbftodet : MS=%d nalpha=%3.2f\n",MS,nalpha);
    //ndets = (int)binom(dNSOMO, nalpha);
    if(NSOMO > 0){
      ndets = (int)binom((double)NSOMO, (double)nalpha);
    }
    else if(NSOMO == 0){
      ndets = 1;
    }
    else printf("Something is wrong in calcMEdetpair\n");

    Tree dettree = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettree.rootNode = malloc(sizeof(Node));
    (*dettree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettree, Isomo, MS, &ndets);

    if(ndets == 1){
    // Initialize transformation matrix
        NBF = 1;
        (*bftodetmatrixptr) = malloc(NBF*ndets*sizeof(double));
        (*rows) = 1;
        (*cols) = 1;

        double *bftodetmatrix = (*bftodetmatrixptr);
        bftodetmatrix[0] = 1.0;

    }
    else{


    int detlist[ndets];
    getDetlistDriver(&dettree, NSOMO, detlist);

    // Prepare BFs
    Tree bftree = (Tree){  .rootNode = NULL, .NBF = -1 };
    bftree.rootNode = malloc(sizeof(Node));
    (*bftree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    generateAllBFs(Isomo, MS, &bftree, &NBF, &NSOMO);

    // Initialize transformation matrix
    //printf("MS=%d NBF=%d ndets=%d NSOMO=%d\n",MS,NBF,ndets,NSOMO);
    assert( NBF > 0);
    assert( ndets > 0);
    (*bftodetmatrixptr) = malloc(NBF*ndets*sizeof(double));
    (*rows) = NBF;
    (*cols) = ndets;

    double *bftodetmatrix = (*bftodetmatrixptr);

    // Build BF to det matrix
    int addI = 0;
    int addJ = 0;
    double rowvec[ndets];
    for(int i=0;i<ndets;i++)
        rowvec[i]=0.0;
    int *BF1 = malloc(MAX_SOMO * sizeof(int));
    int *BF2 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF1 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF2 = malloc(MAX_SOMO * sizeof(int));

    for(int i = 0; i < NBF; i++){
        addI = i;
        getIthBFDriver(&bftree, NSOMO, addI, BF1);
        getBFIndexList(NSOMO, BF1, IdxListBF1);

        // Get ith row
        getbftodetfunction(&dettree, NSOMO, MS, IdxListBF1, rowvec);

        for(int j = 0; j < ndets; j++)
            bftodetmatrix[i*ndets + j] = rowvec[j];

        for(int k=0;k<ndets;k++)
            rowvec[k]=0.0;

        for(int j=0;j<NSOMO;++j){
          BF1[j]=0;
          IdxListBF1[j]=0;
        }
    }

    // Garbage collection
    free(BF1);
    free(IdxListBF1);
    free(BF2);
    free(IdxListBF2);

    }// ndet > 1

}


void convertBFtoDetBasisWithArrayDims(int64_t Isomo, int MS, int rowsmax, int colsmax, int *rows, int *cols, double *bftodetmatrix){

    int NSOMO=0;
    getSetBits(Isomo, &NSOMO);
    int ndets = 0;
    int NBF = 0;
    //double dNSOMO = NSOMO*1.0;
    //double nalpha = (NSOMO + MS)/2.0;
    int nalpha = (NSOMO + MS)/2;
    ndets = (int)binom((double)NSOMO, (double)nalpha);

    Tree dettree = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettree.rootNode = malloc(sizeof(Node));
    (*dettree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettree, Isomo, MS, &ndets);

    //int addr = -1;
    //int inpdet[NSOMO];
    //inpdet[0] = 1;
    //inpdet[1] = 1;
    //inpdet[2] = 1;
    //inpdet[3] = 0;
    //inpdet[4] = 0;
    //inpdet[5] = 0;

    //findAddofDetDriver(&dettree, NSOMO, inpdet, &addr);

    int detlist[ndets];
    getDetlistDriver(&dettree, NSOMO, detlist);

    // Prepare BFs
    Tree bftree = (Tree){  .rootNode = NULL, .NBF = -1 };
    bftree.rootNode = malloc(sizeof(Node));
    (*bftree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    generateAllBFs(Isomo, MS, &bftree, &NBF, &NSOMO);

    // Initialize transformation matrix
    //(*bftodetmatrixptr) = malloc(NBF*ndets*sizeof(double));
    (*rows) = NBF;
    (*cols) = ndets;

    //double *bftodetmatrix = (*bftodetmatrixptr);

    // Build BF to det matrix
    int addI = 0;
    int addJ = 0;
    double rowvec[ndets];
    for(int i=0;i<ndets;i++)
        rowvec[i]=0.0;
    int *BF1 = malloc(MAX_SOMO * sizeof(int));
    int *BF2 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF1 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF2 = malloc(MAX_SOMO * sizeof(int));

    for(int i = 0; i < NBF; i++){
        addI = i;
        getIthBFDriver(&bftree, NSOMO, addI, BF1);
        getBFIndexList(NSOMO, BF1, IdxListBF1);


        // Get ith row
        getbftodetfunction(&dettree, NSOMO, MS, IdxListBF1, rowvec);

        for(int j = 0; j < ndets; j++)
            bftodetmatrix[i*ndets + j] = rowvec[j];

        for(int k=0;k<ndets;k++)
            rowvec[k]=0.0;
    }

    // Garbage collection
    free(BF1);
    free(IdxListBF1);
    free(BF2);
    free(IdxListBF2);

}



void getApqIJMatrixDims(int64_t Isomo, int64_t Jsomo, int64_t MS, int32_t *rowsout, int32_t *colsout){
    int NSOMOI=0;
    int NSOMOJ=0;
    getSetBits(Isomo, &NSOMOI);
    getSetBits(Jsomo, &NSOMOJ);
    int NBFI=0;
    int NBFJ=0;
    getncsfs(NSOMOI, MS, &NBFI);
    getncsfs(NSOMOJ, MS, &NBFJ);
    (*rowsout) = NBFI;
    (*colsout) = NBFJ;
  //exit(0);
}

void getApqIJMatrixDriver(int64_t Isomo, int64_t Jsomo, int orbp, int orbq, int64_t MS, int64_t NMO, double **CSFICSFJApqIJptr, int *rowsout, int *colsout){

    double *overlapMatrixI;
    double *overlapMatrixJ;
    double *orthoMatrixI;
    double *orthoMatrixJ;
    double *bftodetmatrixI;
    double *bftodetmatrixJ;
    double *ApqIJ;
    int NSOMO=0;

    /***********************************
                   Doing I
    ************************************/

    int rowsbftodetI, colsbftodetI;

    convertBFtoDetBasis(Isomo, MS, &bftodetmatrixI, &rowsbftodetI, &colsbftodetI);

    // Fill matrix
    int rowsI = 0;
    int colsI = 0;

    //getOverlapMatrix(Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);
    getOverlapMatrix_withDet(bftodetmatrixI, rowsbftodetI, colsbftodetI, Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);

    orthoMatrixI = malloc(rowsI*colsI*sizeof(double));

    gramSchmidt(overlapMatrixI, rowsI, colsI, orthoMatrixI);

    /***********************************
                   Doing J
    ************************************/


    int rowsbftodetJ, colsbftodetJ;

    convertBFtoDetBasis(Jsomo, MS, &bftodetmatrixJ, &rowsbftodetJ, &colsbftodetJ);

    int rowsJ = 0;
    int colsJ = 0;
    // Fill matrix
    //getOverlapMatrix(Jsomo, MS, &overlapMatrixJ, &rowsJ, &colsJ, &NSOMO);
    getOverlapMatrix_withDet(bftodetmatrixJ, rowsbftodetJ, colsbftodetJ, Jsomo, MS, &overlapMatrixJ, &rowsJ, &colsJ, &NSOMO);

    orthoMatrixJ = malloc(rowsJ*colsJ*sizeof(double));

    gramSchmidt(overlapMatrixJ, rowsJ, colsJ, orthoMatrixJ);

    int rowsA = 0;
    int colsA = 0;

    callcalcMEij(Isomo, Jsomo, orbp, orbq, MS, NMO, &ApqIJ, &rowsA, &colsA);

    // Final ME in BF basis

    // First transform I in bf basis
    double *bfIApqIJ = malloc(rowsbftodetI*colsA*sizeof(double));

    int transA=false;
    int transB=false;
    callBlasMatxMat(bftodetmatrixI, rowsbftodetI, colsbftodetI, ApqIJ, rowsA, colsA, bfIApqIJ, transA, transB);

    // now transform I in csf basis
    double *CSFIApqIJ = malloc(rowsI*colsA*sizeof(double));
    transA = false;
    transB = false;
    callBlasMatxMat(orthoMatrixI, rowsI, colsI, bfIApqIJ, colsI, colsA, CSFIApqIJ, transA, transB);

    // now transform J in BF basis
    double *CSFIbfJApqIJ = malloc(rowsI*rowsbftodetJ*sizeof(double));
    transA = false;
    transB = true;
    callBlasMatxMat(CSFIApqIJ, rowsI, colsA, bftodetmatrixJ, rowsbftodetJ, colsbftodetJ, CSFIbfJApqIJ, transA, transB);

    // now transform J in CSF basis
    (*CSFICSFJApqIJptr) = malloc(rowsI*rowsJ*sizeof(double));
    (*rowsout) = rowsI;
    (*colsout) = rowsJ;

    double *CSFICSFJApqIJ = (*CSFICSFJApqIJptr);
    transA = false;
    transB = true;
    callBlasMatxMat(CSFIbfJApqIJ, rowsI, rowsbftodetJ, orthoMatrixJ, rowsJ, colsJ, CSFICSFJApqIJ, transA, transB);


    // Garbage collection
    free(overlapMatrixI);
    free(overlapMatrixJ);
    free(orthoMatrixI);
    free(orthoMatrixJ);
    free(bftodetmatrixI);
    free(bftodetmatrixJ);
    free(ApqIJ);
    free(bfIApqIJ);
    free(CSFIApqIJ);
    free(CSFIbfJApqIJ);
}

void getApqIJMatrixDriverArrayInp(int64_t Isomo, int64_t Jsomo, int32_t orbp, int32_t orbq, int64_t MS, int64_t NMO, double *CSFICSFJApqIJ, int32_t rowsmax, int32_t colsmax){

    double *overlapMatrixI;
    double *overlapMatrixJ;
    double *orthoMatrixI;
    double *orthoMatrixJ;
    double *bftodetmatrixI;
    double *bftodetmatrixJ;
    double *ApqIJ;
    int NSOMO=0;

    /***********************************
                   Doing I
    ************************************/

    int rowsbftodetI, colsbftodetI;

    //printf(" 1Calling convertBFtoDetBasis Isomo=%ld MS=%ld\n",Isomo,MS);
    convertBFtoDetBasis(Isomo, MS, &bftodetmatrixI, &rowsbftodetI, &colsbftodetI);

    // Fill matrix
    int rowsI = 0;
    int colsI = 0;

    //getOverlapMatrix(Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);
    //printf("Isomo=%ld\n",Isomo);
    getOverlapMatrix_withDet(bftodetmatrixI, rowsbftodetI, colsbftodetI, Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);
    if(Isomo == 0){
      rowsI = 1;
      colsI = 1;
    }

    //printf("Isomo=%ld\n",Isomo);

    orthoMatrixI = malloc(rowsI*colsI*sizeof(double));

    gramSchmidt(overlapMatrixI, rowsI, colsI, orthoMatrixI);

    /***********************************
                   Doing J
    ************************************/

    int rowsbftodetJ, colsbftodetJ;

    //printf(" 2Calling convertBFtoDetBasis Jsomo=%ld MS=%ld\n",Jsomo,MS);
    convertBFtoDetBasis(Jsomo, MS, &bftodetmatrixJ, &rowsbftodetJ, &colsbftodetJ);

    int rowsJ = 0;
    int colsJ = 0;
    // Fill matrix
    //getOverlapMatrix(Jsomo, MS, &overlapMatrixJ, &rowsJ, &colsJ, &NSOMO);
    getOverlapMatrix_withDet(bftodetmatrixJ, rowsbftodetJ, colsbftodetJ, Jsomo, MS, &overlapMatrixJ, &rowsJ, &colsJ, &NSOMO);
    if(Jsomo == 0){
      rowsJ = 1;
      colsJ = 1;
    }

    orthoMatrixJ = malloc(rowsJ*colsJ*sizeof(double));

    gramSchmidt(overlapMatrixJ, rowsJ, colsJ, orthoMatrixJ);

    int rowsA = 0;
    int colsA = 0;

    callcalcMEij(Isomo, Jsomo, orbp, orbq, MS, NMO, &ApqIJ, &rowsA, &colsA);

    // Final ME in BF basis

    // First transform I in bf basis
    double *bfIApqIJ = malloc(rowsbftodetI*colsA*sizeof(double));

    int transA=false;
    int transB=false;
    //printf("1Calling blas\n");
    //printf("rowsA=%d colsA=%d\nrowB=%d colB=%d\n",rowsbftodetI,colsbftodetI,rowsA,colsA);
    callBlasMatxMat(bftodetmatrixI, rowsbftodetI, colsbftodetI, ApqIJ, rowsA, colsA, bfIApqIJ, transA, transB);
    //printf("done\n");

    // now transform I in csf basis
    double *CSFIApqIJ = malloc(rowsI*colsA*sizeof(double));
    transA = false;
    transB = false;
    //printf("2Calling blas\n");
    //printf("rowsA=%d colsA=%d\nrowB=%d colB=%d\n",rowsI,colsI,colsI,colsA);
    callBlasMatxMat(orthoMatrixI, rowsI, colsI, bfIApqIJ, colsI, colsA, CSFIApqIJ, transA, transB);

    // now transform J in BF basis
    double *CSFIbfJApqIJ = malloc(rowsI*rowsbftodetJ*sizeof(double));
    transA = false;
    transB = true;
    //printf("3Calling blas\n");
    //printf("rowsA=%d colsA=%d\nrowB=%d colB=%d\n",rowsI,colsA,rowsbftodetJ,colsbftodetJ);
    callBlasMatxMat(CSFIApqIJ, rowsI, colsA, bftodetmatrixJ, rowsbftodetJ, colsbftodetJ, CSFIbfJApqIJ, transA, transB);

    // now transform J in CSF basis
    //(*CSFICSFJApqIJptr) = malloc(rowsI*rowsJ*sizeof(double));
    //(*rowsout) = rowsI;
    //(*colsout) = rowsJ;

    double *tmpCSFICSFJApqIJ = malloc(rowsI*rowsJ*sizeof(double));
    transA = false;
    transB = true;
    //printf("4Calling blas\n");
    //printf("rowsA=%d colsA=%d\nrowB=%d colB=%d\n",rowsI,rowsbftodetJ,rowsJ,colsJ);
    callBlasMatxMat(CSFIbfJApqIJ, rowsI, rowsbftodetJ, orthoMatrixJ, rowsJ, colsJ, tmpCSFICSFJApqIJ, transA, transB);
    // Transfer to actual buffer in Fortran order
    for(int i = 0; i < rowsI; i++)
        for(int j = 0; j < rowsJ; j++)
            CSFICSFJApqIJ[j*rowsI + i] = tmpCSFICSFJApqIJ[i*rowsJ + j];

    // Garbage collection
    free(overlapMatrixI);
    free(overlapMatrixJ);
    free(orthoMatrixI);
    free(orthoMatrixJ);
    free(bftodetmatrixI);
    free(bftodetmatrixJ);
    free(ApqIJ);
    free(bfIApqIJ);
    free(CSFIApqIJ);
    free(CSFIbfJApqIJ);
    free(tmpCSFICSFJApqIJ);
}

void calculateMETypeSOMOSOMO(int *BF1, int *BF2, int moi, int moj, double *factor, int *phasefactor){

    // Calculate the factor following rules in the table
    // find the type

    if(BF1[moi] == moj){
        // Type I
        (*factor) = sqrt(2.0);
        (*phasefactor) = 1;
    }
    else{
        if(BF1[moi] != moi || BF1[moj] != moj){
            // Type II, III and IV
            (*factor) = 1.0/sqrt(2.0);
            (*phasefactor) =-1;
        }
        else{
            // Type V
            (*factor) = 0.0;
            (*phasefactor) = 1;
        }
    }


}
