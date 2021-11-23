#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bits.h"
#include "parsemap.h"
#include "galois.h"
#define ANF 1
#define pos(i,j) (((1+i)*i)/2 + j)
int *table_red;
typedef struct _mono {
    int w; 
    char ** v;
} mono;

typedef struct _poli {
    int nvars;
    mono ** p;
} poli;


int verify( int ** coef, int * enc, int * result, int npols, int nvars, int w ) {
    int k, i, j, r, n;
    for(k=0; k<npols; k++ ) {
        r = 0;
        // Quadratic terms
        for(i=0;i<nvars;i++) {
            for(j=0;j<=i;j++) {
                    n  = mul( result[i], result[j], w );
                    r = r ^ mul( n, coef[k][pos(i,j)], w );
            }
        }
        // Linear terms
        for( i = 0; i < nvars; i++ ) {
            r = r ^ mul( result[i], coef[k][pos(nvars,i)], w );
        }
        // Constant Terms
        r = r ^ coef[k][pos(nvars,nvars)];
        
        
        if( r !=  enc[k] ) {
            printf("ERROR NOT CONFIRMED POLYNOMIALS!\n" );
            exit( 1 );
        }
    }
    printf("SOLUTION IS CONFIRMED!\n");
    
}


inline char * get( poli * p, int i, int j, int k ) {
    int pos = ((1+i)*i)/2 + j;
    return p->p[pos]->v[k];
} 


mono * new_mono( int w, int idx ) {
    int i;
    mono * r = (mono *) malloc( sizeof(mono) );
    r->w = w;
    r->v = (char **) malloc( sizeof(char *)*( 2*w ) );
    for( i = 0; i < 2*w; i++ ) {
        r->v[i] = (char *) malloc( sizeof(char )*( 4096 ) );
        r->v[i][0] = 0;    
        if( i < w ) {
            #if ANF            
            sprintf(r->v[i], "%d,", idx*w+i+1 );
            #else           
            sprintf(r->v[i], "X%d_%d", idx, i );
            #endif
        }      
    }
    return r;
}

mono * new_mono_clean( int w ) {
    int i;
    mono * r = (mono *) malloc( sizeof(mono) );
    r->w = w;
    r->v = (char **) malloc( sizeof(char *)*( 2*w ) );
    for( i = 0; i < 2*w; i++ ) {
        r->v[i] = (char *) malloc( sizeof(char )*( 4096 ) );
        r->v[i][0] = 0;                  
    }
    return r;
}
mono * mul_mono( mono * m, int v ) {
    int i, j = 0, n;    
    mono * s = new_mono_clean( m->w );
    for(i=0; i < m->w; i++) { 
        j = 0; 
        n = v;       
        while( n != 0 ) {
            if( (n % 2) == 1 ) {
                if( strcmp(s->v[i+j], "" ) == 0 ) 
                    strcat(  s->v[i+j], m->v[i] );               
                else { 
                    #if ANF            
                    sprintf(  s->v[i+j], "%s|%s", s->v[i+j], m->v[i] );               
                    #else
                    sprintf(  s->v[i+j], "%s+%s", s->v[i+j], m->v[i] );
                    #endif
                    // printf("strlen = %d %d\n", strlen( s->v[i+j] ), n );              
                }            
            }
            j++;
            n = n >> 1;
        }
    }
    for(i=0; i< 2*m->w; i++ )
        free( m->v[i] ); 
    free( m );
    return s;
}
mono * new_mono_sq( int w, int idx ) {
    int i;
    mono * r = (mono *) malloc( sizeof(mono) );
    r->w = w;
    r->v = (char **) malloc( sizeof(char *)*(2*w) );
    for( i = 0; i < 2*w; i++ ) {
        r->v[i] = (char *) malloc( sizeof(char )*( 1024 ) );
        r->v[i][0] = 0;    
        strcat(r->v[i], ""); 
    }
    for( i = 0; i < w; i++ ){
        #if ANF   
        sprintf(r->v[2*i], "%d,", idx*w+ i +1);      
        #else     
        sprintf(r->v[2*i], "X%d_%d", idx, i );
        #endif    
    }
    return r;
}

mono * new_mono_mul( int w, int idx_i, int idx_j) {
    int i,j;
    mono * r = (mono *) malloc( sizeof(mono) );
    r->w = w;
    r->v = (char **) malloc( sizeof(char *)*(2*w) );
    for( i = 0; i < 2*w; i++ ) {
        r->v[i] = (char *) malloc( sizeof(char )*( 1024 ) );
        r->v[i][0] = 0;    
        strcat(r->v[i], ""); 
    }
        
    for(i=0;i<w;i++)
        for(j=0;j<w;j++) {
            char t[1024];
            if( strcmp( r->v[i+j], "" ) == 0 ) {
                #if ANF                 
                sprintf(t, "%d,%d,", idx_i*w+ i+1, idx_j*w+ j+1 );                
                #else
                sprintf(t, "X%d_%d*X%d_%d", idx_i, i, idx_j, j );
                #endif
            } else {        
               #if ANF             
               sprintf(t, "|%d,%d,", idx_i*w+ i+1, idx_j*w+ j+1 );   
               #else
               sprintf(t, "+X%d_%d*X%d_%d", idx_i, i, idx_j, j );
               #endif
            }            
            strcat( r->v[i+j], t );             
            //printf("stdlen = %d\n", strlen( r->v[i+j] ) );  
       }
    return r;
} 
 

void reduce( poli * p ){
    int i, j, k, w = p->p[0]->w, nvars = p->nvars, m;  
    for(i=0; i < nvars+1; i++ ) {
        for(j=0; j <= i; j++ ) {
                for(k=w; k <2*w; k++) {
                    // Need reduce
                    if( strcmp( p->p[pos(i,j)]->v[k], "" ) != 0 ) {
                        m = table_red[k];
                        int new_pos = 0;                        
                        while( m != 0 ) {
                            if( (m % 2) == 1 ) {
                                if( strcmp(p->p[pos(i,j)]->v[new_pos], "" ) == 0 )
                                    strcat(p->p[pos(i,j)]->v[new_pos], p->p[pos(i,j)]->v[k]);
                                else  {           
                                    char t[100000];   
                                    #if ANF            
                                    sprintf(t, "%s|%s",p->p[pos(i,j)]->v[new_pos], p->p[pos(i,j)]->v[k]  );
                                    #else 
                                    sprintf(t, "%s+%s",p->p[pos(i,j)]->v[new_pos], p->p[pos(i,j)]->v[k]  );
                                    #endif
                                    strcpy(  p->p[pos(i,j)]->v[new_pos], t );                               
                                } 
                            }
                            new_pos++;
                            m = m >> 1;
                        }
                        p->p[pos(i,j)]->v[k][0] = 0; // clean
                    }   
                }
            
        }    
    }
}

void print_poli( poli * p ) {
    int i, j, k=0;
    for(i=0; i<p->nvars; i++)
        for(j=0; j<=i; j++){
                     
            for( k = 0; k < 2*(p->p[0]->w); k++ ) {
                    if( strcmp( get( p, i, j, k ), "" ) != 0 )
                    printf("(%s)*t^%d + ", get( p, i, j, k ), k );
                            
            }
            printf("\n");
        
        }
            
}

void print_anf( poli * p, int enc, int coef ) {
    FILE * out = fopen("tmp.anf", "a" );
    int i, j, k=0, ff = 1, v = enc, s = coef;
    for( k = 0; k < (p->p[0]->w); k++ ) {
        ff = 1;        
        for(i=0; i<p->nvars+1; i++) {
            for(j=0; j<= i; j++){
                if( strcmp( get( p, i, j, k ), "" ) != 0 ) {
                    if( ff ) {                         
                        fprintf(out, "%s", get( p, i, j, k ) );
                        ff = 0;                        
                    } else {
                        #if ANF                           
                        fprintf(out, "|%s", get( p, i, j, k ) );
                        #else 
                        fprintf(out, "+%s", get( p, i, j, k ) );
                        #endif
                    }    
                }
            }
        }
        #if ANF                           
                          
        if( (v & 1) ^ (s & 1)  )
            fprintf(out, "||\n");
        else
            fprintf(out, "|\n");
        v = v >> 1;
        s = s >> 1;
        #else 
            fprintf(out, "\n");
        #endif
        // printf(" == 0\n");
    }
    fclose( out );         
}

poli * new_poli( int nvars, int w  ) {
    poli * p = (poli *) malloc( sizeof(poli) );
    p->nvars = nvars;
    
    int i,j;
    nvars++;
    p->p = (mono **) malloc( sizeof( mono * ) * ((nvars+1)*nvars)/2);
    nvars--;
    for(i=0; i < nvars; i++ )
        for(j=0; j <=i; j++ )
            if( i == j )                    
                p->p[pos(i,j)] = new_mono_sq( w, i );
            else
                p->p[pos(i,j)] = new_mono_mul( w, i, j );    
    for( i = 0; i < nvars; i++ ) {
        p->p[pos(nvars,i)] = new_mono( w, i );
    }
    p->p[pos(nvars,nvars)] = new_mono_clean( w );
    #if ANF 
    //strcat( p->p[pos(nvars,nvars)]->v[0], "|" ); 
    #else
    strcat( p->p[pos(nvars,nvars)]->v[0], "1+" );
    #endif
    
    reduce( p );
    return p;
} 



/*inline int set( poli * p, int i, int j, int e ) {
    int pos = ((1+i)*i)/2 + j;
    p->p[pos] = e;
}
*/
void create_table( int w ) {
    table_red = (int *) malloc( sizeof( int ) * 2 * w );
    int i, j, a, b;
    for(i=0;i<2*w;i++) {
        table_red[i] = mod( (1 << i), w );
        //printf("%X\n", table_red[i] );
    }
}
int ** read_input( FILE * fin, int * npols, int * nvars, int * w, int ** enc ) {
    int i, j;    
    fscanf( fin, "%d %d %d", npols, nvars, w );
    int nv = *nvars + 1;
    int np = *npols;
    *enc = (int *) malloc( sizeof( int) * np );
    for(i=0; i < np; i++ ) {
        fscanf( fin, "%d", &((*enc)[i]) );
    }
    int ** coefs = (int **) malloc( sizeof( int *) * np );
    for(i=0; i < np; i++ ) {
        coefs[i] =  (int *) malloc( sizeof( int) * ((nv+1)*nv)/2 );    
        for(j=0; j < ((nv+1)*nv)/2; j++)
            fscanf( fin, "%d", &coefs[i][j] );
        
         
    }
    fclose( fin );
    return coefs;
}

int isInMap( int readVar, int nvars, int w, int * map ) {
    int i;
    for( i = 0; i < w*nvars; i++ )
         if( map[i] == readVar )
            return i;
    return -1;
}



int * interpret_result(int nvars, int w, int * map) {
    char f[1024];
    int i, n, res, v = -1, pos, bit, midx;
    int * result = malloc( nvars * sizeof( int ) );
    for(i=0; i < nvars; i++ )
        result[i] = 0;
    
    FILE * r = fopen( "result.sat", "r" );
    fgets( f, 1024, r );
    if( strcmp( f, "SAT\n" ) != 0 ) {
        printf("ERROR UNSAT PROBLEM\n");
        exit( 1 );
    }  else {
        
        while( v != 0 ){
            fscanf( r, "%d", &v );
            midx = isInMap((int)abs(v),nvars,w, map );
            if( midx != -1 ) {
                pos = midx / w;
                bit = midx % w;    
                
                if( v > 0 )
                    setbit( &result[pos], bit );
            }
        }
            
         
    }
    for(i=0; i < nvars; i++ )
        printf("X%d=%d\n", i, result[i] );
    
    return result;
}

void printError() {
   fprintf(stderr, "National Laboratory for Scientific Computing (LNCC) 2014\n");
   fprintf(stderr, "mq2sat [OPTIONS] [FILE]\nOptions are:\n");
   fprintf(stderr, "      -sat sat_solver         Determine a executable to solve boolean equations\n");
   fprintf(stderr, "                              Default is minisat\n");
   fprintf(stderr, "                              Example: $ mq2sat -sat cryptominisat file.mq\n");
   exit( 1 );

}

int main(int argc, char ** argv) {
    int nvars = 2;
    int npols = 128;
    int i;
    int w = 2, *enc;
    system( "rm -rf tmp.anf tmp.cnf");
    if( argc == 1 )
        printError();
    char satExec[1024] = "minisat";
    char fileName[1024];
    fileName[0] = 0;
    for(i=1; i < argc; i++ ) {
        if( strcmp( argv[i], "-sat" ) == 0 ) {
            if( i == argc - 1 )
                printError();
            i++;
            sprintf( satExec ,"%s", argv[i] );
        } else
            strcat( fileName, argv[i] );                        
        
    }
    
    
    FILE * fin = fopen(fileName, "r" );        
    int ** coefs = read_input( fin, &npols, &nvars, &w, &enc );
    int k, j;
    create_table( w ); 
    poli ** p = (poli **) malloc( sizeof( poli *) * npols ); 
    
    for(k=0;k<npols;k++) {
        p[k] = new_poli( nvars, w  );
        fprintf(stderr, "process polynomial .......... %d\n", k);
        for(i=0; i<p[k]->nvars+1; i++)
            for(j=0; j<=i; j++)
                p[k]->p[pos(i,j)] = mul_mono(p[k]->p[pos(i,j)], coefs[k][pos(i,j)] );         
        reduce( p[k] );
        print_anf( p[k], enc[k],  coefs[k][pos(nvars,nvars)] );
    }
    printf("running anf2cnf...\n");
    if( strcmp( satExec, "cryptominisat" ) == 0 )
        system( "anf2cnf --xor --cut=3 tmp.anf > tmp.cnf" );
    else
        system( "anf2cnf --cut=3 tmp.anf > tmp.cnf" );
    
    
    int * map = parse_cnf_map_vars(fopen("tmp.cnf", "r") ,nvars,w );
    printf("running SAT solver...\n");
    
    char runSat[1024];
    sprintf( runSat, "%s tmp.cnf  result.sat", satExec );
    system( runSat );
    int * result = interpret_result(nvars, w, map);
    verify( coefs, enc, result, npols, nvars, w );
    return 0;
}   
