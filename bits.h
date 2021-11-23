#ifndef __BITS_H__
#define __BITS_H__
void resetbit( int * e, int p ) {
   *e &= ~(1 << p); 
}
void setbit( int * e, int p ) {
   *e |= 1 << p; 
}
int bit( int v, int p ) {
    return ((1 << p) & v) ? 1 : 0;
}


int bitlen( int v ) {
    int c = 0;
    while( v != 0 ) {
        v = v >> 1;
        c++;
    } 
    return c;
}
#endif
