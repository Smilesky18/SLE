#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "lu.h"


struct Coordinate {
    int x;
    int y;
    double val;
};
int coordcmp(const void *v1, const void *v2)
{
    struct Coordinate *c1 = (struct Coordinate *) v1;
    struct Coordinate *c2 = (struct Coordinate *) v2;

    if (c1->x != c2->x)
    {
        return (c1->x - c2->x);
    }
    else
    {
        return (c1->y - c2->y);
    }
}

void readMatrix(char *filename, double **val_ptr, int **cols_ptr, 
                int **rowDelimiters_ptr, int *n, int *numRows, int *numCols) 
{
    std::string line;
    char id[128];
    char object[128]; 
    char format[128]; 
    char field[128]; 
    char symmetry[128]; 
    //FILE *mfs = fopen( filename, "r");

    std::ifstream mfs( filename );
    if( !mfs.good() )
    {
        //std::cerr << "Error: unable to open matrix file " << filename << std::endl;
        printf(" Error: unable to open matrix file ");
        exit( 1 );
    }

    int symmetric = 0; 
    int pattern = 0; 
    int field_complex = 0;
    int nRows, nCols, nElements;  
    struct Coordinate *coords;
    // read matrix header
    if( getline( mfs, line ).eof() )
    {
        //std::cerr << "Error: file " << filename << " does not store a matrix" << std::endl;
        printf(" Error: file does not store a matrix ");
        exit( 1 );
    }
    sscanf(line.c_str(), "%s %s %s %s %s", id, object, format, field, symmetry); 
    if (strcmp(object, "matrix") != 0) 
    {
        fprintf(stderr, "Error: file %s does not store a matrix\n", filename); 
        exit(1); 
    }
  
    if (strcmp(format, "coordinate") != 0)
    {
        fprintf(stderr, "Error: matrix representation is dense\n"); 
        exit(1); 
    } 

    if (strcmp(field, "pattern") == 0) 
    {
        pattern = 1; 
    }

	if (strcmp(field, "complex") == 0)
	{
        field_complex = 1;
	}

    if (strcmp(symmetry, "symmetric") == 0) 
    {
        symmetric = 1; 
    }
    while (!getline( mfs, line ).eof() )
    {
        if (line[0] != '%') 
        {
            break; 
        }
    } 
    sscanf(line.c_str(), "%d %d %d", &nRows, &nCols, &nElements); 
    int nElements_padding = (nElements%16 == 0) ? nElements : (nElements + 16)/ 16 * 16;
    int valSize = nElements * sizeof(struct Coordinate);

    if (symmetric) 
    {
        valSize*=2; 
    }                 
    coords = (struct Coordinate*)malloc(valSize);
    int index = 0; 
    float xx99 = 0;
    while (!getline( mfs, line ).eof() )
    {
        if (pattern) 
        {
            sscanf(line.c_str(), "%d %d", &coords[index].x, &coords[index].y); 
            coords[index].val = index%13;
        }
        else if(field_complex)
        {
            // read the value from file
            sscanf(line.c_str(), "%d %d %lf %f", &coords[index].x, &coords[index].y, 
                   &coords[index].val, xx99); 
        }
        else 
        {
            // read the value from file
            sscanf(line.c_str(), "%d %d %lf", &coords[index].x, &coords[index].y, 
                   &coords[index].val); 
        } 
        index++; 
        if (symmetric && coords[index-1].x != coords[index-1].y) 
        {
            coords[index].x = coords[index-1].y; 
            coords[index].y = coords[index-1].x; 
            coords[index].val = coords[index-1].val; 
            index++;
        }

    }  

    nElements = index; 
  
    nElements_padding = (nElements%16 == 0) ? nElements : (nElements + 16)/ 16 * 16;
 
  /*printf("===========================================================================\n");
  printf("=========*********  Informations of the sparse matrix   *********==========\n");
  printf(" Number of Rows is %d\n:", nRows);
  printf(" Number of Columns is %d\n:", nCols);
  printf(" Number of Elements is %d\n:", nElements);
  printf(" After Alignment : %d\n", nElements_padding);
  printf("===========================================================================\n");
  printf("............ Converting the Raw matrix to CSR .................\n");*/ 
  /*for (int qq = index; qq <nElements_padding; qq ++)
  {
            coords[qq].x = coords[index - 1].x;
            coords[qq].y = coords[index - 1].y;
            coords[qq].val = 0;      
  }*/
  for ( int i = 0; i < nElements; i++ )
  {
    coords[i].x = coords[i].x - 1;
    coords[i].y = coords[i].y - 1;
  }
  qsort(coords, nElements, sizeof(struct Coordinate), coordcmp); 
  /*for ( int i = 0; i < nElements; i++ )
  {
    printf(" coords[%d].x = %d\n", i, coords[i].x);
    printf(" coords[%d].y = %d\n", i, coords[i].y);
    printf(" coords[%d].val = %f\n", i, coords[i].val);
  }*/


    // create CSR data structures
    *n = nElements; 
    *numRows = nRows; 
    *numCols = nCols; 
//    *val_ptr = new floatType[nElements_padding]; 
//    *cols_ptr = new int[nElements_padding];
//    *rowDelimiters_ptr = new int[nRows+2]; 

    *val_ptr = (double *)malloc(sizeof(double) * nElements); 
    *cols_ptr = (int *)malloc(sizeof(int) * nElements);
    *rowDelimiters_ptr = (int *)malloc(sizeof(int) * (nRows+2));

    double *val = *val_ptr; 
    int *cols = *cols_ptr; 
    int *rowDelimiters = *rowDelimiters_ptr; 
    rowDelimiters[0] = 0; 
    int r=0; 
    int i=0;
    for (i=0; i<nElements; i++) 
    {
        while (coords[i].x != r) 
        {
            rowDelimiters[++r] = i; 
        }
        val[i] = coords[i].val; 
        cols[i] = coords[i].y;    
    }
    for(int k = r + 1; k<=(nRows+1); k++)
    {
       rowDelimiters[k] = i; 
    }
    r = 0; 
    free(coords); 
}