/* Agent X -- meet-in-the-middle over low-Hamming-weight subsets of the 256 ladder points.
   Field: p = 2^256 - 2^32 - 977 (secp256k1 prime).  Curve y^2 = x^3 + b, a = 0.
   Modes:  table  -> emit low-64 bits of x( sum_{i in A} 2^i G ) for all |A| in [1..SMAX]
           scan   -> for all |B| = SZ, compute T - sum_{i in B} 2^i G, look its key up
           bitmap -> build the 2^32-bit prefilter from a sorted key file
           find   -> locate which subset produced a given key (hit post-processing)
           selftest
*/
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef u64 fe[4];

static const u64 CC = 0x1000003D1ULL;                 /* p = 2^256 - CC */
static const u64 PP[4] = {0xFFFFFFFEFFFFFC2FULL, ~0ULL, ~0ULL, ~0ULL};

static inline void fe_copy(u64*r,const u64*a){r[0]=a[0];r[1]=a[1];r[2]=a[2];r[3]=a[3];}
static inline int  fe_iszero(const u64*a){return (a[0]|a[1]|a[2]|a[3])==0;}
static inline int  fe_eq(const u64*a,const u64*b){return a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2]&&a[3]==b[3];}
static inline int  fe_ge_p(const u64*a){
    if(a[3]!=PP[3]) return a[3]>PP[3];
    if(a[2]!=PP[2]) return a[2]>PP[2];
    if(a[1]!=PP[1]) return a[1]>PP[1];
    return a[0]>=PP[0];
}
static inline void fe_subp(u64*a){ /* a -= p, assumes a>=p */
    u128 b=0; for(int i=0;i<4;i++){u128 v=(u128)a[i]-PP[i]-b; a[i]=(u64)v; b=(v>>64)&1;}
}
static inline void fe_add(u64*r,const u64*a,const u64*b){
    u64 c=0,t[4];
    for(int i=0;i<4;i++){u128 v=(u128)a[i]+b[i]+c; t[i]=(u64)v; c=(u64)(v>>64);}
    if(c){ /* 2^256 == CC mod p */
        u128 v=(u128)t[0]+CC; t[0]=(u64)v; u64 c2=(u64)(v>>64);
        for(int i=1;i<4&&c2;i++){u128 w=(u128)t[i]+c2; t[i]=(u64)w; c2=(u64)(w>>64);}
    }
    if(fe_ge_p(t)) fe_subp(t);
    fe_copy(r,t);
}
static inline void fe_sub(u64*r,const u64*a,const u64*b){
    u64 t[4]; u128 br=0;
    for(int i=0;i<4;i++){u128 v=(u128)a[i]-b[i]-br; t[i]=(u64)v; br=(v>>64)&1;}
    if(br){ u64 c=0; for(int i=0;i<4;i++){u128 v=(u128)t[i]+PP[i]+c; t[i]=(u64)v; c=(u64)(v>>64);} }
    fe_copy(r,t);
}
static inline void fe_neg(u64*r,const u64*a){ if(fe_iszero(a)){r[0]=r[1]=r[2]=r[3]=0;return;} fe_sub(r,PP,a); }

static inline void fe_mul(u64*r,const u64*a,const u64*b){
    u64 t[8]; t[0]=t[1]=t[2]=t[3]=t[4]=t[5]=t[6]=t[7]=0;
    for(int i=0;i<4;i++){
        u64 carry=0;
        for(int j=0;j<4;j++){
            u128 v=(u128)a[i]*b[j]+t[i+j]+carry;
            t[i+j]=(u64)v; carry=(u64)(v>>64);
        }
        t[i+4]=carry;
    }
    /* fold: res = t[0..3] + t[4..7]*CC */
    u64 res[4]; u64 carry=0;
    for(int j=0;j<4;j++){ u128 v=(u128)t[4+j]*CC + t[j] + carry; res[j]=(u64)v; carry=(u64)(v>>64); }
    /* carry < 2^34 ; fold again */
    while(carry){
        u128 v=(u128)carry*CC + res[0]; res[0]=(u64)v; u64 c2=(u64)(v>>64); carry=0;
        for(int i=1;i<4&&c2;i++){u128 w=(u128)res[i]+c2; res[i]=(u64)w; c2=(u64)(w>>64);}
        carry=c2;
    }
    if(fe_ge_p(res)) fe_subp(res);
    fe_copy(r,res);
}
static void fe_inv(u64*r,const u64*a){
    static const u64 e[4]={0xFFFFFFFEFFFFFC2DULL,~0ULL,~0ULL,~0ULL};  /* p-2 */
    u64 res[4]={1,0,0,0}, base[4]; fe_copy(base,a); int started=0;
    for(int i=255;i>=0;i--){
        if(started) fe_mul(res,res,res);
        if((e[i>>6]>>(i&63))&1){ if(!started){fe_copy(res,base);started=1;} else fe_mul(res,res,base); }
    }
    fe_copy(r,res);
}


/* --- checked read-only mmap.  A missing/empty table must ABORT with a clear message, never
   segfault on first dereference: an unchecked MAP_FAILED once cost another agent six scans. --- */
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
static const void* xmap_ro(const char*path, size_t*nbytes){
    int fd = open(path, O_RDONLY);
    if(fd < 0){ fprintf(stderr,"FATAL: cannot open '%s' (missing?)\n", path); exit(2); }
    struct stat st;
    if(fstat(fd,&st) != 0 || st.st_size == 0){
        fprintf(stderr,"FATAL: '%s' is empty or unstattable\n", path); exit(2); }
    void*m = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if(m == MAP_FAILED){ fprintf(stderr,"FATAL: mmap('%s', %lld bytes) failed\n",
                                 path,(long long)st.st_size); exit(2); }
    if(nbytes) *nbytes = (size_t)st.st_size;
    return m;
}
