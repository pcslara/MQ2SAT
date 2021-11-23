#ifndef __GALOIS_H__
#define __GALOIS_H__

int prim_poly[33] = 
{ 0, 
/*  1 */     2, 
/*  2 */    07,
/*  3 */    013,
/*  4 */    023,
/*  5 */    045,
/*  6 */    0103,
/*  7 */    0211,
/*  8 */    0435,
/*  9 */    01021,
/* 10 */    02011,
/* 11 */    04005,
/* 12 */    010123,
/* 13 */    020033,
/* 14 */    042103,
/* 15 */    0100003,
/* 16 */    0210013,
/* 17 */    0400011,
/* 18 */    01000201,
/* 19 */    02000047,
/* 20 */    04000011,
/* 21 */    010000005,
/* 22 */    020000003,
/* 23 */    040000041,
/* 24 */    0100000207,
/* 25 */    0200000011,
/* 26 */    0400000107,
/* 27 */    01000000047,
/* 28 */    02000000011,
/* 29 */    04000000005,
/* 30 */    010040000007,
/* 31 */    020000000011, 
/* 32 */    00020000007 };  /* Really 40020000007, but we're omitting the high order bit */


int mul( int a, int b, int w ) {
    int i,j,r=0, d;    
    for(i=0;i<w;i++) {
        for(j=0;j<w;j++) {
            d = bit(r, i + j ) ^ (bit( a, i )&bit( b, j ));
            if( d )
                setbit( &r, i + j );
            else
                resetbit( &r, i + j );                 
        }
    }
    return mod( r, w );
}

int mod( int a, int w ) {
    int ret = a;
    int modular = prim_poly[w]; 
    while( 1 ) {
        int diff = bitlen( ret ) - bitlen( prim_poly[w] ); 
        if( diff < 0 )
            return ret;
        modular = prim_poly[w];
        modular = modular << diff;
        ret = ret ^ modular;
    }
}
#endif
