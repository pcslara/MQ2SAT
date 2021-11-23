#ifndef __PARSEMAP_H__ 
#define __PARSEMAP_H__

void get_values( int * idx, int *var, char * in ) {
    int iIn = 0;
    int iOut = 0;
    int i;
    char strNum[20];    
    while( isdigit( in[iIn] ) == 0 )
        iIn++;

    i = 0;
    while( isdigit( in[iIn] ) ) {
        strNum[i] = in[iIn];
        i++; iIn++;      
    }
    strNum[i] = 0;
    *var = atoi( strNum ); 
    //printf("FIND VAR := %d\n", *var ); 
    while( isdigit( in[iIn] ) == 0 )
        iIn++;
    i = 0;
    while( isdigit( in[iIn] ) ) {
        strNum[i] = in[iIn];
        i++; iIn++;      
    }
    strNum[i] = 0;
    *idx = atoi( strNum ) - 1;    
    //printf("FIND IDX := %d\n", *idx );
}

int * parse_cnf_map_vars( FILE * fCNF, int nvars, int w ) {
    char line[1024];
    int * map = (int *) malloc( sizeof( int ) * nvars * w);
    int idx, var, n;
    for(n=0; n < nvars*w; n++ )
        map[n] = -1;
        
    while( fgets( line, 1024, fCNF ) != NULL ) {
        if( line[0] == 'c' &&  strstr( line, "x" ) != NULL && strstr( line, "*" ) == NULL ) {
            get_values( &idx, &var, line );
            map[idx] = var;
        }
    }
    return map;

}

#endif
