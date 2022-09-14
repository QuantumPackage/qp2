#include <assert.h>
#include "tree_utils.h"

void buildTree(Tree *bftree,
               Node **inode,
               int isomo,
               int izeros,
               int icpl,
               int NSOMOMax,
               int MSmax){

    // Find the maximum parallel couplings 0
    //      the maximum anti-parallel couplings 1
    int zeromax = MSmax + (NSOMOMax-MSmax)/2;
    int onemax = NSOMOMax - zeromax;

    // Exit condition
    if(isomo > NSOMOMax || icpl < 0 || izeros > zeromax ) return;

    // If we find a valid BF assign its address
    if(isomo == NSOMOMax && icpl == MSmax){
        (*inode)->addr = bftree->rootNode->addr;
        bftree->rootNode->addr += 1;
        return;
    }

    // Call 0 branch
    if(izeros+1 <= zeromax){
      if(((*inode)->C0) == NULL){
          ((*inode)->C0) = malloc(sizeof(Node));
          (*(*inode)->C0) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = *inode, .addr = -1, .cpl = 0, .iSOMO = isomo };
          buildTree(bftree, &(*inode)->C0, isomo+1, izeros+1, icpl+1, NSOMOMax, MSmax);
      }
      else buildTree(bftree, &(*inode)->C0, isomo+1, izeros+1, icpl+1, NSOMOMax, MSmax);
    }

    // Call 1 branch
    if(icpl-1 >=0){
      if(((*inode)->C1) == NULL){
          ((*inode)->C1) = malloc(sizeof(Node));
          (*(*inode)->C1) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = *inode, .addr = -1, .cpl = 1, .iSOMO = isomo };
          buildTree(bftree, &(*inode)->C1, isomo+1, izeros+0, icpl-1, NSOMOMax, MSmax);
      }
      else buildTree(bftree, &(*inode)->C1, isomo+1, izeros+0, icpl-1, NSOMOMax, MSmax);
    }

    return;
}

void buildTreeDriver(Tree *bftree, int NSOMO, int MS, int *NBF){
    int isomo = 0; // counts the total number of SOMO's
    int izeros= 0; // Counts the total number of parallel coupings (i.e. 0's)
    int icpl  = 0; // keep track of the ith ms (cannot be -ve)
    int addr  = 0; // Counts the total BF's

    assert(bftree->rootNode->addr == 0);
    buildTree(bftree, &(bftree->rootNode), isomo, izeros, icpl, NSOMO, MS);

    *NBF = bftree->rootNode->addr;
}


void getIthBF(Node *inode, int isomo, bool foundBF, int NSOMOMax, int getaddr, int *vecBF){
    // Exit condition
    if(foundBF) return;
    if(isomo > NSOMOMax) return;
    if(inode == NULL) return;

    if(isomo == NSOMOMax){
        if(inode->addr == getaddr){
            int i;
            for(i = NSOMOMax-1; i > -1; i--){
                vecBF[i] = inode->cpl;
                inode = inode->PREV;
            }
            foundBF = true;
            return;
        }
    }

    // Recurse to C0
    if(inode->C0 != NULL){
        getIthBF(inode->C0, isomo+1, foundBF, NSOMOMax, getaddr, vecBF);
    }
    // Recurse to C1
    if(inode->C1 != NULL){
        getIthBF(inode->C1, isomo+1, foundBF, NSOMOMax, getaddr, vecBF);
    }

    return;
}

void getIthBFDriver(Tree *bftree, int NSOMOMax, int getaddr, int *vecBF){
    int isomo = 0;
    bool foundBF = false;
    getIthBF((bftree->rootNode), isomo, foundBF, NSOMOMax, getaddr, vecBF);
}

void genDets(Tree *dettree,
               Node **inode,
               int isomo,
               int izeros,
               int icpl,
               int NSOMOMax,
               int MSmax){

    // Find the maximum parallel couplings 0
    //      the maximum anti-parallel couplings 1
    int zeromax = MSmax + (NSOMOMax-MSmax)/2;
    int onemax = NSOMOMax - zeromax;

    // Exit condition
    if(isomo > NSOMOMax || izeros > zeromax || abs(icpl) > onemax) return;

    // If we find a valid BF assign its address
    if(isomo == NSOMOMax){
        (*inode)->addr = dettree->rootNode->addr;
        dettree->rootNode->addr += 1;
        return;
    }

    // Call 0 branch
    if(((*inode)->C0) == NULL && izeros+1 <= zeromax){
        ((*inode)->C0) = malloc(sizeof(Node));
        (*(*inode)->C0) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = *inode, .addr = -1, .cpl = 0, .iSOMO = isomo };
        genDets(dettree, &(*inode)->C0, isomo+1, izeros+1, icpl+0, NSOMOMax, MSmax);
    }
    else genDets(dettree, &(*inode)->C0, isomo+1, izeros+1, icpl+0, NSOMOMax, MSmax);

    // Call 1 branch
    if(((*inode)->C1) == NULL && abs(icpl+1) <= onemax){
        ((*inode)->C1) = malloc(sizeof(Node));
        (*(*inode)->C1) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = *inode, .addr = -1, .cpl = 1, .iSOMO = isomo };
        genDets(dettree, &(*inode)->C1, isomo+1, izeros+0, icpl+1, NSOMOMax, MSmax);
    }
    else genDets(dettree, &(*inode)->C1, isomo+1, izeros+0, icpl+1, NSOMOMax, MSmax);

    return;
}

void genDetsDriver(Tree *dettree, int NSOMO, int MS, int *Ndets){
    int isomo = 0; // counts the total number of SOMO's
    int izeros= 0; // Counts the total number of parallel coupings (i.e. 0's)
    int icpl  = 0; // keep track of the ith ms (cannot be -ve)
    int addr  = 0; // Counts the total BF's

    genDets(dettree, &(dettree->rootNode), isomo, izeros, icpl, NSOMO, MS);

    *Ndets = dettree->rootNode->addr;
}

void getIthDet(Node *inode, int isomo, bool foundBF, int NSOMOMax, int getaddr, int *vecBF){
    // Exit condition
    if(foundBF) return;
    if(isomo > NSOMOMax) return;
    if(inode == NULL) return;

    if(isomo == NSOMOMax){
        if(inode->addr == getaddr){
            int i;
            for(i = NSOMOMax-1; i > -1; i--){
                vecBF[i] = inode->cpl;
                inode = inode->PREV;
            }
            foundBF = true;
            return;
        }
    }

    // Recurse to C0
    if(inode->C0 != NULL){
        getIthDet(inode->C0, isomo+1, foundBF, NSOMOMax, getaddr, vecBF);
    }
    // Recurse to C1
    if(inode->C1 != NULL){
        getIthDet(inode->C1, isomo+1, foundBF, NSOMOMax, getaddr, vecBF);
    }

    return;
}

void getIthDetDriver(Tree *dettree, int NSOMOMax, int getaddr, int *vecBF){
    int isomo = 0;
    bool foundBF = false;
    getIthDet((dettree->rootNode), isomo, foundBF, NSOMOMax, getaddr, vecBF);
}

void findAddofDet(Node *inode, int isomo, bool foundDet, int NSOMOMax, int *inpdet, int *addr){
    // Exit condition
    if(foundDet) return;
    if(isomo == NSOMOMax){
        foundDet = true;
        *addr = inode->addr;
        return;
    }

    // Recurse to C0
    if(inpdet[isomo] == 0){
        if(inode->C0 != NULL){
            findAddofDet(inode->C0, isomo+1, foundDet, NSOMOMax, inpdet, addr);
        }
        else{
            *addr = -1;
            return;
        }
    }
    else{
        // Recurse to C1
        if(inode->C1 != NULL){
            findAddofDet(inode->C1, isomo+1, foundDet, NSOMOMax, inpdet, addr);
        }
        else{
            *addr = -1;
            return;
        }
    }

    return;
}

void findAddofDetDriver(Tree *dettree, int NSOMOMax, int *inpdet, int *addr){
    *addr = -1;
    int isomo = 0;
    bool foundDet = false;
    findAddofDet((dettree->rootNode), isomo, foundDet, NSOMOMax, inpdet, addr);
}

void getDetlist(Node *inode, int isomo, int NSOMOMax, int *vecBF, int *detlist){
    // Exit condition
    if(isomo > NSOMOMax) return;
    if(inode == NULL) return;

    if(isomo == NSOMOMax){
        int idet=0;
        int k;
        for(k=0;k<NSOMOMax;k++){
            if(vecBF[k] == 1) idet = idet | (1<<(NSOMOMax-1-k));
        }
        detlist[inode->addr]=idet;
        return;
    }


    // Recurse to C0
    if(inode->C0 != NULL){
        vecBF[isomo] = 0;
        getDetlist(inode->C0, isomo+1, NSOMOMax, vecBF, detlist);
    }
    // Recurse to C1
    if(inode->C1 != NULL){
        vecBF[isomo] = 1;
        getDetlist(inode->C1, isomo+1, NSOMOMax, vecBF, detlist);
    }

    return;
}

void getDetlistDriver(Tree *dettree, int NSOMOMax, int *detlist){
    int isomo = 0;
    int vecBF[NSOMOMax];
    if(dettree->rootNode->addr > 1) getDetlist((dettree->rootNode), isomo, NSOMOMax, vecBF, detlist);
}

void genDetBasis(Tree *dettree, int Isomo, int MS, int *ndets){

    int NSOMO=0;
    getSetBits(Isomo, &NSOMO);
    genDetsDriver(dettree, NSOMO, MS, ndets);
    // Closed shell case
    if(NSOMO==0) (*ndets) = 1;

}

void callBlasMatxMat(double *A, int rowA, int colA, double *B, int rowB, int colB, double *C, bool transA, bool transB){
    int m = rowA;
    int k = colA;
    int n = colB;
    double alpha = 1.0;
    double beta  = 0.0;
    int val = 0;

    if (transA) val |= 0x1;
    if (transB) val |= 0x2;

    switch (val) {
        case 0: // notransA, notransB
            m = rowA;
            n = colB;
            k = colA;
            f_dgemm('N', 'N', n, m, k, alpha, B, n, A, k, beta, C, n);
            break;
        case 1: // transA, notransB
            m = colA;
            n = colB;
            k = rowA;
            f_dgemm('N', 'T', n, m, k, alpha, B, n, A, colA, beta, C, n);
            break;
        case 2: // notransA, transB
            //m = rowA;
            //n = rowB;
            //k = colB;
            m = rowA;
            n = rowB;
            k = colA;
            f_dgemm('T', 'N', n, m, k, alpha, B, colB, A, k, beta, C, n);
            break;
        case 3: // transA, transB
            m = colA;
            n = rowB;
            k = rowA;
            f_dgemm('T', 'T', n, m, k, alpha, B, colB, A, colA, beta, C, n);
            break;
        default:
            printf("Impossible !!!!\n");
            break;
    }
}

void printRealMatrix(double *orthoMatrix, int rows, int cols){
  int i,j;
  for(i=0;i<rows;++i){
    for(j=0;j<cols;++j){
      printf(" %3.5f ",orthoMatrix[i*cols + j]);
    }
    printf("\n");
  }
}
