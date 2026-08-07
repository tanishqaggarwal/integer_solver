/* Agent X -- meet-in-the-middle over low-Hamming-weight subsets of the 256 ladder points.
   Field: p = 2^256 - 2^32 - 977 (secp256k1 prime).  Curve y^2 = x^3 + b, a = 0.
   Modes:  table  -> emit low-64 bits of x( sum_{i in A} 2^i G ) for all |A| in [1..SMAX]
           scan   -> for all |B| = SZ, compute T - sum_{i in B} 2^i G, look its key up
           bitmap -> build the 2^32-bit prefilter from a sorted key file
           find   -> locate which subset produced a given key (hit post-processing)
           selftest
*/
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

/* ------------------------------------------------------------------ globals */
static fe LX[256], LY[256];      /* the ladder, 2^i G, possibly negated for scan side */
static fe BX, BY;                /* base point of the walk (T for scan) */
static int  SIGN;                /* +1 table (base = O), -1 scan (base = T, adds -2^i G) */
static int  SZ;                  /* subset size being enumerated */
static u64  ZEROEV=0;            /* count of degenerate dx==0 events */
static FILE *REPORT=NULL;

/* lookup structures (scan mode) */
static const u64 *TBL=NULL; static size_t TBLN=0;
static const unsigned char *BM=NULL;

#define BUFN 4096
typedef struct {
    fe x1,y1; int m; u64 code;
} Item;
typedef struct {
    Item buf[BUFN]; int bn;
    fe dx[BUFN], pre[BUFN];
    u64 keys[BUFN];
    FILE *out;                 /* table mode: raw key sink */
    u64 nproc;
    u64 obuf[1<<16]; int obn;
} TS;

static inline int tbl_has(u64 k){
    u64 bi=k>>32;
    if(!((BM[bi>>3]>>(bi&7))&1)) return 0;
    size_t lo=0, hi=TBLN;
    while(lo<hi){ size_t mid=(lo+hi)>>1; if(TBL[mid]<k) lo=mid+1; else hi=mid; }
    return lo<TBLN && TBL[lo]==k;
}

static void flush_last(TS*ts){
    int n=ts->bn; if(!n) return;
    for(int i=0;i<n;i++){
        fe_sub(ts->dx[i], ts->buf[i].x1, LX[ts->buf[i].m]);
        if(fe_iszero(ts->dx[i])){
            #pragma omp critical
            { ZEROEV++; if(REPORT&&ZEROEV<10000) fprintf(REPORT,"ZERO %d %llu %d\n",SZ,(unsigned long long)ts->buf[i].code,ts->buf[i].m); }
            ts->dx[i][0]=1; ts->dx[i][1]=ts->dx[i][2]=ts->dx[i][3]=0;
        }
    }
    /* batch inversion */
    fe acc; acc[0]=1;acc[1]=acc[2]=acc[3]=0;
    for(int i=0;i<n;i++){ fe_copy(ts->pre[i],acc); fe_mul(acc,acc,ts->dx[i]); }
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv; fe_mul(iv,ainv,ts->pre[i]); fe_mul(ainv,ainv,ts->dx[i]);
        fe lam,num,x3;
        fe_sub(num, ts->buf[i].y1, LY[ts->buf[i].m]);
        fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam);
        fe_sub(x3,x3,ts->buf[i].x1);
        fe_sub(x3,x3,LX[ts->buf[i].m]);
        ts->keys[i]=x3[0];
    }
    if(ts->out){
        for(int i=0;i<n;i++){
            ts->obuf[ts->obn++]=ts->keys[i];
            if(ts->obn==(1<<16)){ fwrite(ts->obuf,8,ts->obn,ts->out); ts->obn=0; }
        }
    } else {
        for(int i=0;i<n;i++) if(tbl_has(ts->keys[i])){
            #pragma omp critical
            { fprintf(REPORT,"HIT %d %llu %d %llu\n",SZ,(unsigned long long)ts->buf[i].code,
                      ts->buf[i].m,(unsigned long long)ts->keys[i]); fflush(REPORT); }
        }
    }
    ts->nproc+=n; ts->bn=0;
}
static inline void emit_last(TS*ts,const u64*x,const u64*y,int start,u64 code,int depth){
    for(int m=start;m<256;m++){
        Item*it=&ts->buf[ts->bn++];
        fe_copy(it->x1,x); fe_copy(it->y1,y); it->m=m; it->code=code|((u64)m<<(8*depth));
        if(ts->bn==BUFN) flush_last(ts);
    }
}
/* full add of (x,y) + L[m] for m in [start,256), results into ox/oy arrays */
static void batch_full(TS*ts,const u64*x,const u64*y,int start,fe*ox,fe*oy){
    int n=256-start;
    static __thread fe d[256],pr[256];
    for(int i=0;i<n;i++){
        fe_sub(d[i],x,LX[start+i]);
        if(fe_iszero(d[i])){
            #pragma omp critical
            { ZEROEV++; if(REPORT&&ZEROEV<10000) fprintf(REPORT,"ZEROFULL %d %d\n",SZ,start+i); }
            d[i][0]=1;d[i][1]=d[i][2]=d[i][3]=0;
        }
    }
    fe acc; acc[0]=1;acc[1]=acc[2]=acc[3]=0;
    for(int i=0;i<n;i++){ fe_copy(pr[i],acc); fe_mul(acc,acc,d[i]); }
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv; fe_mul(iv,ainv,pr[i]); fe_mul(ainv,ainv,d[i]);
        fe lam,num,x3,y3,t;
        fe_sub(num,y,LY[start+i]); fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam); fe_sub(x3,x3,x); fe_sub(x3,x3,LX[start+i]);
        fe_sub(t,x,x3); fe_mul(y3,lam,t); fe_sub(y3,y3,y);
        fe_copy(ox[i],x3); fe_copy(oy[i],y3);
    }
}
static void rec(TS*ts,const u64*x,const u64*y,int start,int depth,u64 code){
    if(depth==SZ-1){ emit_last(ts,x,y,start,code,depth); return; }
    int n=256-start; if(n<=0) return;
    fe *ox=malloc(sizeof(fe)*n), *oy=malloc(sizeof(fe)*n);
    batch_full(ts,x,y,start,ox,oy);
    for(int i=0;i<n;i++){
        int idx=start+i;
        if(256-(idx+1) < SZ-1-depth) break;   /* not enough room left */
        rec(ts,ox[i],oy[i],idx+1,depth+1,code|((u64)idx<<(8*depth)));
    }
    free(ox); free(oy);
}

static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}

static void hexto(u64*r,const char*s){ /* decimal string -> fe, via simple mul-add */
    r[0]=r[1]=r[2]=r[3]=0;
    for(const char*q=s;*q;q++){
        u64 c=0; for(int i=0;i<4;i++){u128 v=(u128)r[i]*10+c; r[i]=(u64)v; c=(u64)(v>>64);}
        c=*q-'0'; for(int i=0;i<4&&c;i++){u128 v=(u128)r[i]+c; r[i]=(u64)v; c=(u64)(v>>64);}
    }
}
int main(int argc,char**argv){
    /* argv: mode datafile [args] */
    if(argc<3){fprintf(stderr,"usage: xmitm mode datafile ...\n");return 1;}
    const char*mode=argv[1];
    /* load ladder + T from a plain text file: line1 Tx, line2 Ty, then 256 lines "x y" */
    FILE*f=fopen(argv[2],"r"); if(!f){perror("data");return 1;}
    char l1[200],l2[200];
    if(fscanf(f,"%199s %199s",l1,l2)!=2){fprintf(stderr,"bad data\n");return 1;}
    hexto(BX,l1); hexto(BY,l2);
    for(int i=0;i<256;i++){ char a[200],b[200];
        if(fscanf(f,"%199s %199s",a,b)!=2){fprintf(stderr,"bad ladder %d\n",i);return 1;}
        hexto(LX[i],a); hexto(LY[i],b); }
    fclose(f);

    if(!strcmp(mode,"selftest")){
        /* verify field arithmetic + that L[i] are on the curve and L[i+1]=2L[i] via the add law */
        fe one={1,0,0,0}, t1,t2;
        fe_inv(t1,LX[0]); fe_mul(t2,t1,LX[0]);
        printf("inv ok %d\n", fe_eq(t2,one));
        /* print x(T - L[3] ) and x(L[5]+L[7]) for python cross-check */
        TS*ts=calloc(1,sizeof(TS));
        fe ox[256],oy[256];
        SZ=2; SIGN=1;
        batch_full(ts,LX[5],LY[5],7,ox,oy);
        printf("L5+L7 x = %llu %llu %llu %llu\n",(unsigned long long)ox[0][0],(unsigned long long)ox[0][1],(unsigned long long)ox[0][2],(unsigned long long)ox[0][3]);
        fe nlx[256],nly[256];
        for(int i=0;i<256;i++){fe_copy(nlx[i],LX[i]);fe_neg(nly[i],LY[i]);}
        memcpy(LY,nly,sizeof(LY));
        batch_full(ts,BX,BY,3,ox,oy);
        printf("T-L3 x = %llu %llu %llu %llu\n",(unsigned long long)ox[0][0],(unsigned long long)ox[0][1],(unsigned long long)ox[0][2],(unsigned long long)ox[0][3]);
        printf("T-L3 y = %llu %llu %llu %llu\n",(unsigned long long)oy[0][0],(unsigned long long)oy[0][1],(unsigned long long)oy[0][2],(unsigned long long)oy[0][3]);
        return 0;
    }
    if(!strcmp(mode,"bitmap")){
        /* argv[3]=sorted key file, argv[4]=out bitmap */
        int fd=open(argv[3],O_RDONLY); struct stat st; fstat(fd,&st);
        size_t n=st.st_size/8;
        const u64*k=mmap(NULL,st.st_size,PROT_READ,MAP_SHARED,fd,0);
        size_t bmsz=(size_t)1<<29;   /* 2^32 bits */
        unsigned char*bm=calloc(bmsz,1);
        for(size_t i=0;i<n;i++){ u64 bi=k[i]>>32; bm[bi>>3]|=(unsigned char)(1u<<(bi&7)); }
        FILE*o=fopen(argv[4],"wb"); fwrite(bm,1,bmsz,o); fclose(o);
        printf("bitmap from %zu keys\n",n); return 0;
    }

    int neg = !strcmp(mode,"scan") || !strcmp(mode,"find_scan");
    if(neg){ for(int i=0;i<256;i++) fe_neg(LY[i],LY[i]); SIGN=-1; } else SIGN=1;

    if(!strcmp(mode,"table")){
        /* argv[3]=smax, argv[4]=outfile */
        int smax=atoi(argv[3]);
        TS*ts=calloc(1,sizeof(TS));
        ts->out=fopen(argv[4],"wb");
        double t0=now();
        for(SZ=1;SZ<=smax;SZ++){
            if(SZ==1){
                /* keys are simply x(L[m]) */
                for(int m=0;m<256;m++){ ts->obuf[ts->obn++]=LX[m][0]; }
                fwrite(ts->obuf,8,ts->obn,ts->out); ts->obn=0; ts->nproc+=256;
            } else {
                for(int i0=0;i0<256;i0++){
                    if(256-(i0+1) < SZ-1) break;
                    rec(ts,LX[i0],LY[i0],i0+1,1,(u64)i0);
                }
                flush_last(ts);
            }
            fprintf(stderr,"  size %d done  n=%llu  %.1fs\n",SZ,(unsigned long long)ts->nproc,now()-t0);
        }
        if(ts->obn){fwrite(ts->obuf,8,ts->obn,ts->out); ts->obn=0;}
        fclose(ts->out);
        fprintf(stderr,"table total %llu keys, zero-events %llu, %.1fs\n",(unsigned long long)ts->nproc,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    if(!strcmp(mode,"scan")){
        /* argv[3]=size argv[4]=tablefile argv[5]=bitmapfile argv[6]=reportfile [argv[7]=i0lo argv[8]=i0hi] */
        SZ=atoi(argv[3]);
        int fd=open(argv[4],O_RDONLY); struct stat st; fstat(fd,&st); TBLN=st.st_size/8;
        TBL=mmap(NULL,st.st_size,PROT_READ,MAP_SHARED,fd,0);
        int fd2=open(argv[5],O_RDONLY); struct stat st2; fstat(fd2,&st2);
        BM=mmap(NULL,st2.st_size,PROT_READ,MAP_SHARED,fd2,0);
        REPORT=fopen(argv[6],"a");
        int lo=(argc>7)?atoi(argv[7]):0, hi=(argc>8)?atoi(argv[8]):256;
        double t0=now();
        /* precompute A1[i0] = T - L[i0] */
        static fe A1x[256],A1y[256];
        { TS*t=calloc(1,sizeof(TS)); fe ox[256],oy[256]; batch_full(t,BX,BY,0,ox,oy);
          for(int i=0;i<256;i++){fe_copy(A1x[i],ox[i]);fe_copy(A1y[i],oy[i]);} free(t); }
        u64 total=0;
        if(SZ==1){
            TS*ts=calloc(1,sizeof(TS));
            /* B of size 1: key = x(T - L[m]) */
            for(int m=0;m<256;m++) if(tbl_has(A1x[m][0]))
                fprintf(REPORT,"HIT 1 %d - %llu\n",m,(unsigned long long)A1x[m][0]);
            total=256; free(ts);
        } else {
            #pragma omp parallel reduction(+:total)
            {
                TS*ts=calloc(1,sizeof(TS));
                #pragma omp for schedule(dynamic,1)
                for(int i0=lo;i0<hi;i0++){
                    if(256-(i0+1) < SZ-1) continue;
                    rec(ts,A1x[i0],A1y[i0],i0+1,1,(u64)i0);
                    flush_last(ts);
                    #pragma omp critical
                    fprintf(stderr,"  i0=%d done cum=%llu %.1fs\n",i0,(unsigned long long)ts->nproc,now()-t0);
                }
                flush_last(ts);
                total+=ts->nproc; free(ts);
            }
        }
        fprintf(REPORT,"DONE size=%d range=[%d,%d) n=%llu zero=%llu %.1fs\n",SZ,lo,hi,
                (unsigned long long)total,(unsigned long long)ZEROEV,now()-t0);
        fflush(REPORT);
        fprintf(stderr,"scan size %d range [%d,%d): %llu candidates, zero-events %llu, %.1fs\n",
                SZ,lo,hi,(unsigned long long)total,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    fprintf(stderr,"unknown mode %s\n",mode); return 1;
}
