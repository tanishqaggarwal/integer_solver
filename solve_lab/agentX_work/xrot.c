/* Agent X -- |S| = 10 via the 128-rotation splitting system.
   A_j = {(j+t) mod 256 : t=0..127},  B_j = complement.  Every 10-set has some j in [0,128)
   with |S n A_j| = |S n B_j| = 5 (discrete continuity: f(j+1)-f(j) in {-1,0,1} and
   f(j)+f(j+128) = 10).  So: table = 5-subsets of A_j, scan = 5-subsets of B_j.
   modes: build <data> <j> <lo> <hi> <out>
          merge <out> <in1> ...
          bitmap <sorted> <out>
          scan  <data> <j> <sortedtbl> <bm> <rep> <lo> <hi>                                  */
#include "xfield.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define NP 128
#define SZ 5
static fe LX[NP],LY[NP],BX,BY,GX[256],GY[256];
static u64 ZEROEV=0; static FILE*REPORT=NULL;
static const u64*TBL=NULL; static size_t TBLN=0; static const unsigned char*BM=NULL;
#define BUFN 4096
typedef struct { fe x1,y1; int m; u64 code; } Item;
typedef struct { Item buf[BUFN]; int bn; fe dx[BUFN],pre[BUFN]; u64 keys[BUFN];
                 FILE*out; u64 nproc; u64 obuf[1<<16]; int obn; } TS;
static inline int tbl_has(u64 k){
    u64 bi=k>>32; if(!((BM[bi>>3]>>(bi&7))&1)) return 0;
    size_t lo=0,hi=TBLN; while(lo<hi){size_t mid=(lo+hi)>>1; if(TBL[mid]<k)lo=mid+1; else hi=mid;}
    return lo<TBLN && TBL[lo]==k;
}
static void flush_last(TS*ts){
    int n=ts->bn; if(!n) return;
    for(int i=0;i<n;i++){ fe_sub(ts->dx[i],ts->buf[i].x1,LX[ts->buf[i].m]);
        if(fe_iszero(ts->dx[i])){
            #pragma omp critical
            { ZEROEV++; if(REPORT&&ZEROEV<10000) fprintf(REPORT,"ZERO %llu %d\n",(unsigned long long)ts->buf[i].code,ts->buf[i].m); }
            ts->dx[i][0]=1;ts->dx[i][1]=ts->dx[i][2]=ts->dx[i][3]=0; } }
    fe acc={1,0,0,0};
    for(int i=0;i<n;i++){fe_copy(ts->pre[i],acc);fe_mul(acc,acc,ts->dx[i]);}
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv;fe_mul(iv,ainv,ts->pre[i]);fe_mul(ainv,ainv,ts->dx[i]);
        fe lam,num,x3;
        fe_sub(num,ts->buf[i].y1,LY[ts->buf[i].m]); fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam); fe_sub(x3,x3,ts->buf[i].x1); fe_sub(x3,x3,LX[ts->buf[i].m]);
        ts->keys[i]=x3[0];
    }
    if(ts->out){ for(int i=0;i<n;i++){ ts->obuf[ts->obn++]=ts->keys[i];
            if(ts->obn==(1<<16)){fwrite(ts->obuf,8,ts->obn,ts->out);ts->obn=0;} } }
    else { for(int i=0;i<n;i++) if(tbl_has(ts->keys[i])){
            #pragma omp critical
            { fprintf(REPORT,"HIT %llu %d %llu\n",(unsigned long long)ts->buf[i].code,
                      ts->buf[i].m,(unsigned long long)ts->keys[i]); fflush(REPORT);} } }
    ts->nproc+=n; ts->bn=0;
}
static inline void emit_last(TS*ts,const u64*x,const u64*y,int start,u64 code,int depth){
    for(int m=start;m<NP;m++){ Item*it=&ts->buf[ts->bn++];
        fe_copy(it->x1,x); fe_copy(it->y1,y); it->m=m; it->code=code|((u64)m<<(8*depth));
        if(ts->bn==BUFN) flush_last(ts); }
}
static void batch_full(TS*ts,const u64*x,const u64*y,int start,fe*ox,fe*oy){
    int n=NP-start; static __thread fe d[NP],pr[NP];
    for(int i=0;i<n;i++){ fe_sub(d[i],x,LX[start+i]);
        if(fe_iszero(d[i])){
            #pragma omp critical
            { ZEROEV++; if(REPORT&&ZEROEV<10000) fprintf(REPORT,"ZEROFULL %d\n",start+i); }
            d[i][0]=1;d[i][1]=d[i][2]=d[i][3]=0; } }
    fe acc={1,0,0,0};
    for(int i=0;i<n;i++){fe_copy(pr[i],acc);fe_mul(acc,acc,d[i]);}
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv;fe_mul(iv,ainv,pr[i]);fe_mul(ainv,ainv,d[i]);
        fe lam,num,x3,y3,t;
        fe_sub(num,y,LY[start+i]); fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam); fe_sub(x3,x3,x); fe_sub(x3,x3,LX[start+i]);
        fe_sub(t,x,x3); fe_mul(y3,lam,t); fe_sub(y3,y3,y);
        fe_copy(ox[i],x3); fe_copy(oy[i],y3);
    }
}
static void rec(TS*ts,const u64*x,const u64*y,int start,int depth,u64 code){
    if(depth==SZ-1){ emit_last(ts,x,y,start,code,depth); return; }
    int n=NP-start; if(n<=0) return;
    fe *ox=malloc(sizeof(fe)*n),*oy=malloc(sizeof(fe)*n);
    batch_full(ts,x,y,start,ox,oy);
    for(int i=0;i<n;i++){ int idx=start+i;
        if(NP-(idx+1) < SZ-1-depth) break;
        rec(ts,ox[i],oy[i],idx+1,depth+1,code|((u64)idx<<(8*depth))); }
    free(ox);free(oy);
}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void hexto(u64*r,const char*s){r[0]=r[1]=r[2]=r[3]=0;
    for(const char*q=s;*q;q++){u64 c=0;for(int i=0;i<4;i++){u128 v=(u128)r[i]*10+c;r[i]=(u64)v;c=(u64)(v>>64);}
        c=*q-'0';for(int i=0;i<4&&c;i++){u128 v=(u128)r[i]+c;r[i]=(u64)v;c=(u64)(v>>64);} } }
static void loadladder(const char*fn){
    FILE*f=fopen(fn,"r"); char a[200],b[200];
    if(fscanf(f,"%199s %199s",a,b)!=2) exit(1); hexto(BX,a); hexto(BY,b);
    for(int i=0;i<256;i++){ if(fscanf(f,"%199s %199s",a,b)!=2) exit(1); hexto(GX[i],a); hexto(GY[i],b); }
    fclose(f);
}
int main(int argc,char**argv){
    const char*mode=argv[1];
    if(!strcmp(mode,"merge")){
        int nf=argc-3; FILE**f=malloc(sizeof(FILE*)*nf); u64*cur=malloc(8*nf); int*live=malloc(4*nf);
        for(int i=0;i<nf;i++){ f[i]=fopen(argv[3+i],"rb"); live[i]=fread(&cur[i],8,1,f[i])==1; }
        FILE*o=fopen(argv[2],"wb"); u64*ob=malloc(8*(1<<16)); int on=0; u64 n=0;
        for(;;){ int best=-1;
            for(int i=0;i<nf;i++) if(live[i] && (best<0 || cur[i]<cur[best])) best=i;
            if(best<0) break;
            ob[on++]=cur[best]; n++;
            if(on==(1<<16)){fwrite(ob,8,on,o);on=0;}
            live[best]=fread(&cur[best],8,1,f[best])==1;
        }
        if(on) fwrite(ob,8,on,o);
        fclose(o); fprintf(stderr,"merged %llu keys\n",(unsigned long long)n); return 0;
    }
    if(!strcmp(mode,"bitmap")){
        size_t nb; const u64*k=xmap_ro(argv[2],&nb); size_t n=nb/8;
        size_t bmsz=(size_t)1<<29; unsigned char*bm=calloc(bmsz,1);
        for(size_t i=0;i<n;i++){u64 bi=k[i]>>32; bm[bi>>3]|=(unsigned char)(1u<<(bi&7));}
        FILE*o=fopen(argv[3],"wb"); fwrite(bm,1,bmsz,o); fclose(o);
        fprintf(stderr,"bitmap from %zu keys\n",n); return 0;
    }
    loadladder(argv[2]);
    int j=atoi(argv[3]);
    if(!strcmp(mode,"build")){
        for(int t=0;t<NP;t++){ int pos=(j+t)&255; fe_copy(LX[t],GX[pos]); fe_copy(LY[t],GY[pos]); }
        int lo=atoi(argv[4]),hi=atoi(argv[5]);
        TS*ts=calloc(1,sizeof(TS)); ts->out=fopen(argv[6],"wb");
        double t0=now();
        for(int i0=lo;i0<hi;i0++){ if(NP-(i0+1)<SZ-1) break; rec(ts,LX[i0],LY[i0],i0+1,1,(u64)i0); }
        flush_last(ts);
        if(ts->obn){fwrite(ts->obuf,8,ts->obn,ts->out);ts->obn=0;}
        fclose(ts->out);
        fprintf(stderr,"build j=%d [%d,%d) n=%llu zero=%llu %.1fs\n",j,lo,hi,
                (unsigned long long)ts->nproc,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    if(!strcmp(mode,"scan")){
        for(int t=0;t<NP;t++){ int pos=(j+128+t)&255; fe_copy(LX[t],GX[pos]); fe_neg(LY[t],GY[pos]); }
        size_t nb; TBL=xmap_ro(argv[4],&nb); TBLN=nb/8;
        size_t nb2; BM=xmap_ro(argv[5],&nb2);
        if(nb2 != ((size_t)1<<29)){fprintf(stderr,"FATAL: bitmap '%s' is %zu bytes, expected %zu\n",argv[5],nb2,(size_t)1<<29);exit(2);}
        REPORT=fopen(argv[6],"a");
        int lo=atoi(argv[7]),hi=atoi(argv[8]);
        double t0=now();
        static fe A1x[NP],A1y[NP];
        { TS*t=calloc(1,sizeof(TS)); fe ox[NP],oy[NP]; batch_full(t,BX,BY,0,ox,oy);
          for(int i=0;i<NP;i++){fe_copy(A1x[i],ox[i]);fe_copy(A1y[i],oy[i]);} free(t); }
        TS*ts=calloc(1,sizeof(TS));
        for(int i0=lo;i0<hi;i0++){ if(NP-(i0+1)<SZ-1) break; rec(ts,A1x[i0],A1y[i0],i0+1,1,(u64)i0); }
        flush_last(ts);
        fprintf(REPORT,"DONE rot=%d range=[%d,%d) n=%llu zero=%llu %.1fs\n",j,lo,hi,
                (unsigned long long)ts->nproc,(unsigned long long)ZEROEV,now()-t0);
        fflush(REPORT);
        fprintf(stderr,"scan j=%d [%d,%d) n=%llu zero=%llu %.1fs\n",j,lo,hi,
                (unsigned long long)ts->nproc,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    return 1;
}
