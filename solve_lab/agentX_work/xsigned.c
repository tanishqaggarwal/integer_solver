/* Agent X -- SIGNED-DIGIT meet-in-the-middle:  k = sum_j eps_j * 2^{e_j},  eps_j in {+1,-1}.
   Strictly contains the unsigned weight class, and also contains low RUN-LENGTH k
   (a run of ones from bit b to bit a is 2^{a+1} - 2^b -- two signed terms).

   Signed ladder index s in [0,512):  position m = s>>1, sign = (s&1) ? -1 : +1.
   Exponents strictly increase, so after using s the next index must be >= ((s>>1)+1)<<1.

   modes: table <data> <smax> <out>            a-term table, LEADING SIGN FIXED POSITIVE
                                               (WLOG: x(P) = x(-P) absorbs global negation)
          scan  <data> <sz> <tbl> <bm> <rep> [s0lo s0hi]
*/
#include "xfield.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static fe SX[512],SY[512],BX,BY;
static int SZ; static u64 ZEROEV=0; static FILE*REPORT=NULL;
static const u64*TBL=NULL; static size_t TBLN=0; static const unsigned char*BM=NULL;
#define BUFN 4096
typedef struct { fe x1,y1; int s; u64 code; } Item;
typedef struct { Item buf[BUFN]; int bn; fe dx[BUFN],pre[BUFN]; u64 keys[BUFN];
                 FILE*out; u64 nproc; u64 obuf[1<<16]; int obn; } TS;
static inline int nxt(int s){ return ((s>>1)+1)<<1; }
static inline int tbl_has(u64 k){
    u64 bi=k>>32; if(!((BM[bi>>3]>>(bi&7))&1)) return 0;
    size_t lo=0,hi=TBLN; while(lo<hi){size_t mid=(lo+hi)>>1; if(TBL[mid]<k)lo=mid+1; else hi=mid;}
    return lo<TBLN && TBL[lo]==k;
}
static void flush_last(TS*ts){
    int n=ts->bn; if(!n) return;
    for(int i=0;i<n;i++){
        fe_sub(ts->dx[i],ts->buf[i].x1,SX[ts->buf[i].s]);
        if(fe_iszero(ts->dx[i])){
            #pragma omp critical
            { ZEROEV++; if(REPORT&&ZEROEV<10000) fprintf(REPORT,"ZERO %d %llu %d\n",SZ,(unsigned long long)ts->buf[i].code,ts->buf[i].s); }
            ts->dx[i][0]=1;ts->dx[i][1]=ts->dx[i][2]=ts->dx[i][3]=0;
        }
    }
    fe acc={1,0,0,0};
    for(int i=0;i<n;i++){fe_copy(ts->pre[i],acc);fe_mul(acc,acc,ts->dx[i]);}
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv;fe_mul(iv,ainv,ts->pre[i]);fe_mul(ainv,ainv,ts->dx[i]);
        fe lam,num,x3;
        fe_sub(num,ts->buf[i].y1,SY[ts->buf[i].s]); fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam); fe_sub(x3,x3,ts->buf[i].x1); fe_sub(x3,x3,SX[ts->buf[i].s]);
        ts->keys[i]=x3[0];
    }
    if(ts->out){ for(int i=0;i<n;i++){ ts->obuf[ts->obn++]=ts->keys[i];
            if(ts->obn==(1<<16)){fwrite(ts->obuf,8,ts->obn,ts->out);ts->obn=0;} } }
    else { for(int i=0;i<n;i++) if(tbl_has(ts->keys[i])){
            #pragma omp critical
            { fprintf(REPORT,"HIT %d %llu %d %llu\n",SZ,(unsigned long long)ts->buf[i].code,
                      ts->buf[i].s,(unsigned long long)ts->keys[i]); fflush(REPORT); } } }
    ts->nproc+=n; ts->bn=0;
}
static inline void emit_last(TS*ts,const u64*x,const u64*y,int start,u64 code,int depth){
    for(int s=start;s<512;s++){
        Item*it=&ts->buf[ts->bn++];
        fe_copy(it->x1,x); fe_copy(it->y1,y); it->s=s; it->code=code|((u64)s<<(16*depth));
        if(ts->bn==BUFN) flush_last(ts);
    }
}
static void batch_full(TS*ts,const u64*x,const u64*y,int start,fe*ox,fe*oy){
    int n=512-start; static __thread fe d[512],pr[512];
    for(int i=0;i<n;i++){ fe_sub(d[i],x,SX[start+i]);
        if(fe_iszero(d[i])){
            #pragma omp critical
            { ZEROEV++; if(REPORT&&ZEROEV<10000) fprintf(REPORT,"ZEROFULL %d %d\n",SZ,start+i); }
            d[i][0]=1;d[i][1]=d[i][2]=d[i][3]=0; } }
    fe acc={1,0,0,0};
    for(int i=0;i<n;i++){fe_copy(pr[i],acc);fe_mul(acc,acc,d[i]);}
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv;fe_mul(iv,ainv,pr[i]);fe_mul(ainv,ainv,d[i]);
        fe lam,num,x3,y3,t;
        fe_sub(num,y,SY[start+i]); fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam); fe_sub(x3,x3,x); fe_sub(x3,x3,SX[start+i]);
        fe_sub(t,x,x3); fe_mul(y3,lam,t); fe_sub(y3,y3,y);
        fe_copy(ox[i],x3); fe_copy(oy[i],y3);
    }
}
static void rec(TS*ts,const u64*x,const u64*y,int start,int depth,u64 code){
    if(depth==SZ-1){ emit_last(ts,x,y,start,code,depth); return; }
    int n=512-start; if(n<=0) return;
    fe *ox=malloc(sizeof(fe)*n),*oy=malloc(sizeof(fe)*n);
    batch_full(ts,x,y,start,ox,oy);
    for(int i=0;i<n;i++){ int s=start+i;
        if(nxt(s) > 512-2*(SZ-1-depth)) break;
        rec(ts,ox[i],oy[i],nxt(s),depth+1,code|((u64)s<<(16*depth)));
    }
    free(ox);free(oy);
}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void hexto(u64*r,const char*s){r[0]=r[1]=r[2]=r[3]=0;
    for(const char*q=s;*q;q++){u64 c=0;for(int i=0;i<4;i++){u128 v=(u128)r[i]*10+c;r[i]=(u64)v;c=(u64)(v>>64);}
        c=*q-'0';for(int i=0;i<4&&c;i++){u128 v=(u128)r[i]+c;r[i]=(u64)v;c=(u64)(v>>64);} } }
int main(int argc,char**argv){
    const char*mode=argv[1];
    FILE*f=fopen(argv[2],"r"); char a[200],b[200];
    if(fscanf(f,"%199s %199s",a,b)!=2)return 1; hexto(BX,a); hexto(BY,b);
    for(int i=0;i<256;i++){ if(fscanf(f,"%199s %199s",a,b)!=2)return 1;
        hexto(SX[2*i],a); hexto(SY[2*i],b);
        fe_copy(SX[2*i+1],SX[2*i]); fe_neg(SY[2*i+1],SY[2*i]); }
    fclose(f);
    double t0=now();
    if(!strcmp(mode,"table")){
        int smax=atoi(argv[3]); TS*ts=calloc(1,sizeof(TS)); ts->out=fopen(argv[4],"wb");
        for(SZ=1;SZ<=smax;SZ++){
            /* leading index EVEN only (sign of the lowest exponent fixed to +1) */
            for(int s0=0;s0<512;s0+=2){
                if(SZ==1){ ts->obuf[ts->obn++]=SX[s0][0]; ts->nproc++;
                    if(ts->obn==(1<<16)){fwrite(ts->obuf,8,ts->obn,ts->out);ts->obn=0;} continue; }
                if(nxt(s0) > 512-2*(SZ-1)) break;
                rec(ts,SX[s0],SY[s0],nxt(s0),1,(u64)s0);
            }
            flush_last(ts);
            fprintf(stderr,"  a=%d done n=%llu %.1fs\n",SZ,(unsigned long long)ts->nproc,now()-t0);
        }
        if(ts->obn){fwrite(ts->obuf,8,ts->obn,ts->out);ts->obn=0;}
        fclose(ts->out);
        fprintf(stderr,"signed table total %llu keys zero=%llu %.1fs\n",(unsigned long long)ts->nproc,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    if(!strcmp(mode,"scan")){
        SZ=atoi(argv[3]);
        size_t nb; TBL=xmap_ro(argv[4],&nb); TBLN=nb/8;
        size_t nb2; BM=xmap_ro(argv[5],&nb2);
        if(nb2 != ((size_t)1<<29)){fprintf(stderr,"FATAL: bitmap '%s' is %zu bytes, expected %zu\n",argv[5],nb2,(size_t)1<<29);exit(2);}
        REPORT=fopen(argv[6],"a");
        int lo=(argc>7)?atoi(argv[7]):0, hi=(argc>8)?atoi(argv[8]):512;
        static fe A1x[512],A1y[512];
        { TS*t=calloc(1,sizeof(TS)); fe ox[512],oy[512]; batch_full(t,BX,BY,0,ox,oy);
          for(int i=0;i<512;i++){fe_copy(A1x[i],ox[i]);fe_copy(A1y[i],oy[i]);} free(t); }
        u64 total=0;
        if(SZ==1){ TS*ts=calloc(1,sizeof(TS));
            for(int s=lo;s<hi;s++) if(tbl_has(A1x[s][0]))
                fprintf(REPORT,"HIT 1 %d - %llu\n",s,(unsigned long long)A1x[s][0]);
            total=hi-lo; free(ts);
        } else {
            #pragma omp parallel reduction(+:total)
            { TS*ts=calloc(1,sizeof(TS));
              #pragma omp for schedule(dynamic,1)
              for(int s0=lo;s0<hi;s0++){
                  if(nxt(s0) > 512-2*(SZ-1)) continue;
                  rec(ts,A1x[s0],A1y[s0],nxt(s0),1,(u64)s0);
                  flush_last(ts);
                  #pragma omp critical
                  fprintf(stderr,"  s0=%d done cum=%llu %.1fs\n",s0,(unsigned long long)ts->nproc,now()-t0);
              }
              flush_last(ts); total+=ts->nproc; free(ts); }
        }
        fprintf(REPORT,"DONE signed sz=%d range=[%d,%d) n=%llu zero=%llu %.1fs\n",SZ,lo,hi,
                (unsigned long long)total,(unsigned long long)ZEROEV,now()-t0);
        fflush(REPORT);
        fprintf(stderr,"signed scan sz=%d [%d,%d): %llu cand zero=%llu %.1fs\n",SZ,lo,hi,
                (unsigned long long)total,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    return 1;
}
