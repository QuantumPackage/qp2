#ifndef TREE_UTILS_H
#define TREE_UTILS_H
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

typedef struct bin_node Node;
typedef struct bin_tree Tree;
struct bin_node {
    Node *C0;
    Node *C1;
    Node *PREV;
    int addr;
    int cpl;
    int iSOMO;
};

struct bin_tree {
    Node *rootNode;
    int NBF;
};

void f_dgemm(const char transb, const char transa, const int n, const int m, const int k,
             const double alpha, const double* B, const int ldb, const double* A,
             const int lda, const double beta, double* C, const int ldc);



#define MAX_SOMO 32

void buildTreeDriver(Tree *bftree, int NSOMO, int MS, int *NBF);

void buildTree(Tree *bftree, Node **inode, int isomo, int izeros, int icpl, int NSOMOMax, int MSmax);

void printTreeDriver(Tree *bftree, int NSOMOMax);
void printTree(Node *bftree, int isomo, int NSOMOMax, int *vecBF);

void getIthBF(Node *node, int isomo, bool foundBF, int NSOMOMax, int getaddr, int *vecBF);
void getIthBFDriver(Tree *bftree, int NSOMOMax, int getaddr, int *vecBF);

void getBFIndexList(int NSOMO, int *BF1, int *IdxListBF1);
void getIslands(int NSOMO, int *BF1, int *BF2, int *nislands, int *phasefactor);

void generateAllBFs(int64_t Isomo, int64_t MS, Tree *bftree, int *NBF, int *NSOMO);
void getSetBits(int64_t n, int *nsetbits);
void getOverlapMatrix(int64_t Isomo, int64_t MS, double **overlapMatrixptr, int *rows, int *cols, int *NSOMOout);
void getOverlapMatrix_withDet(double *bftodetmatrixI, int rowsbftodetI, int colsbftodetI, int64_t Isomo, int64_t MS, double **overlapMatrixI, int *rowsI, int *colsI, int *NSOMO);
void gramSchmidt_qp(double *overlapMatrix, int rows, int cols, double *orthoMatrix);
void gramSchmidt(double *overlapMatrix, int rows, int cols, double *orthoMatrix);


void calculateMETypeSOMOSOMO(int *BF1, int *BF2, int moi, int moj, double *factor, int *phasefactor);
void getOneElMETypeSOMOSOMO(int64_t Isomo, int64_t Jsomos, int moi, int moj, int MS, double **oneElMatrixElementsptr, int *rows, int *cols);

/***********************

Determinant Tree utils
***********************/


void genDets(Tree *dettree,
               Node **inode,
               int isomo,
               int izeros,
               int icpl,
               int NSOMOMax,
               int MSmax);
void genDetsDriver(Tree *dettree, int NSOMO, int MS, int *Ndets);

void getIthDet(Node *inode, int isomo, bool foundBF, int NSOMOMax, int getaddr, int *vecBF);
void getIthDetDriver(Tree *dettree, int NSOMOMax, int getaddr, int *vecBF);
void getDetlistDriver(Tree *dettree, int NSOMOMax, int *detlist);
void findAddofDet(Node *inode, int isomo, bool foundDet, int NSOMOMax, int *inpdet, int *addr);
void findAddofDetDriver(Tree *dettree, int NSOMOMax, int *inpdet, int *addr);


/************************/

void genDetBasis(Tree *dettree, int Isomo, int MS, int *ndets);
void getbftodetfunction(Tree *dettree, int NSOMO, int MS, int *BF1, double *rowvec);
void convertBFtoDetBasis(int64_t Isomo, int MS, double **bftodetmatrixptr, int *rows, int *cols);

// Misc utils
void int_to_bin_digit(int64_t in, int count, int* out);
void printRealMatrix(double *orthoMatrix, int rows, int cols);
void callBlasMatxMat(double *A, int rowA, int colA, double *B, int rowB, int colB, double *C, bool transA, bool transB);
void get_phase_cfg_to_qp_inpList(int *inpdet, int NSOMO, int *phaseout);
void get_phase_cfg_to_qp_inpInt(int inpdet, double *phaseout);

#endif
