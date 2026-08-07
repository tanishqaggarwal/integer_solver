/* Agent AA -- OFFSET-SHIFTED signed-digit meet-in-the-middle.
 *
 *   k*G = T   <=>   (k-c)*G = T - c*G
 *
 * Derived from agent X's xsigned.c (same field code, same signed-ladder indexing,
 * same batched-inversion inner loop).  Two changes:
 *
 *   1. The base point B read from the data file is now T - c*G for an offset c.
 *      Nothing else in the algorithm changes -- that is the whole point of the angle.
 *   2. The table is SHARDED by the top 3 bits of the key into 8 files, so an
 *      a<=4 signed table (1,409,460,736 keys / 11.3 GB) can be built and sorted
 *      with ~1.5 GB of RAM instead of 11.3 GB.  Because shard order agrees with key
 *      order, per-shard binary search is exact.
 *
 * Signed ladder index s in [0,512): exponent = s>>1, sign = (s&1) ? -1 : +1.
 *
 * modes:
 *   table  <data> <smax> <prefix> [nthreads]   -> <prefix>.<shard>.<thread>, leading sign +
 *   bitmap <prefix> <bm.bin>
 *   scan   <data> <sz> <prefix> <bm> <rep> [s0lo s0hi]
 */
#include "aa_field.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define NSH 8
#define SHIFT_SH 61
#define OBN (1<<16)

static fe SX[512],SY[512],BX,BY;
static int SZ; static u64 ZEROEV=0; static FILE*REPORT=NULL;
static const u64*TBL[NSH]; static size_t TBLN[NSH]; static const unsigned char*BM=NULL;
#define BUFN 4096
typedef struct { fe x1,y1; int s; u64 code; } Item;
typedef struct { Item buf[BUFN]; int bn; fe dx[BUFN],pre[BUFN]; u64 keys[BUFN];
                 FILE*out[NSH]; u64 obuf[NSH][OBN]; int obn[NSH];
                 u64 nproc; } TS;

static inline int nxt(int s){ return ((s>>1)+1)<<1; }

static inline int tbl_has(u64 k){
    u64 bi=k>>32; if(!((BM[bi>>3]>>(bi&7))&1)) return 0;
    int sh=(int)(k>>SHIFT_SH);
    const u64*t=TBL[sh]; size_t lo=0,hi=TBLN[sh];
    while(lo<hi){size_t mid=(lo+hi)>>1; if(t[mid]<k)lo=mid+1; else hi=mid;}
    return lo<TBLN[sh] && t[lo]==k;
}
static inline void put(TS*ts,u64 k){
    int sh=(int)(k>>SHIFT_SH);
    ts->obuf[sh][ts->obn[sh]++]=k;
    if(ts->obn[sh]==OBN){ fwrite(ts->obuf[sh],8,OBN,ts->out[sh]); ts->obn[sh]=0; }
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
    if(ts->out[0]){ for(int i=0;i<n;i++) put(ts,ts->keys[i]); }
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

static void load_data(const char*path){
    FILE*f=fopen(path,"r"); if(!f){fprintf(stderr,"no data %s\n",path);exit(2);}
    char a[200],b[200];
    if(fscanf(f,"%199s %199s",a,b)!=2){fprintf(stderr,"bad data\n");exit(2);} hexto(BX,a); hexto(BY,b);
    for(int i=0;i<256;i++){ if(fscanf(f,"%199s %199s",a,b)!=2){fprintf(stderr,"bad data\n");exit(2);}
        hexto(SX[2*i],a); hexto(SY[2*i],b);
        fe_copy(SX[2*i+1],SX[2*i]); fe_neg(SY[2*i+1],SY[2*i]); }
    fclose(f);
}

int main(int argc,char**argv){
    if(argc<2) return 1;
    const char*mode=argv[1];
    double t0=now();

    if(!strcmp(mode,"bitmap")){
        const char*pref=argv[2];
        unsigned char*bm=calloc(1ull<<29,1);      /* 2^32 bits */
        u64 tot=0;
        for(int sh=0;sh<NSH;sh++){
            char fn[512]; snprintf(fn,sizeof fn,"%s.%d",pref,sh);
            FILE*f=fopen(fn,"rb"); if(!f){fprintf(stderr,"missing %s\n",fn);return 2;}
            static u64 buf[1<<16]; size_t n;
            while((n=fread(buf,8,1<<16,f))>0){ for(size_t i=0;i<n;i++){ u64 bi=buf[i]>>32; bm[bi>>3]|=1u<<(bi&7);} tot+=n; }
            fclose(f);
        }
        FILE*o=fopen(argv[3],"wb"); fwrite(bm,1,1ull<<29,o); fclose(o);
        u64 set=0; for(u64 i=0;i<(1ull<<29);i++){unsigned char v=bm[i]; while(v){set+=v&1;v>>=1;}}
        fprintf(stderr,"bitmap: %llu keys, %llu of 2^32 bits set (%.2f%%) %.1fs\n",
                (unsigned long long)tot,(unsigned long long)set,100.0*set/4294967296.0,now()-t0);
        return 0;
    }

    load_data(argv[2]);

    if(!strcmp(mode,"table")){
        int smax=atoi(argv[3]); const char*pref=argv[4];
        int nt=(argc>5)?atoi(argv[5]):1;
        u64 grand=0;
        #pragma omp parallel num_threads(nt) reduction(+:grand)
        {
          int tid=0;
          #ifdef _OPENMP
          tid=omp_get_thread_num();
          #endif
          TS*ts=calloc(1,sizeof(TS));
          for(int sh=0;sh<NSH;sh++){ char fn[512];
              snprintf(fn,sizeof fn,"%s.%d.%d",pref,sh,tid); ts->out[sh]=fopen(fn,"wb"); }
          for(int sz=1;sz<=smax;sz++){
              #pragma omp barrier
              SZ=sz;
              #pragma omp for schedule(dynamic,1)
              for(int s0=0;s0<512;s0+=2){       /* leading index EVEN: sign of lowest exponent = + */
                  if(sz==1){ put(ts,SX[s0][0]); ts->nproc++; continue; }
                  if(nxt(s0) > 512-2*(sz-1)) continue;
                  rec(ts,SX[s0],SY[s0],nxt(s0),1,(u64)s0);
                  flush_last(ts);
              }
              flush_last(ts);
              #pragma omp critical
              fprintf(stderr,"  a=%d tid=%d cum=%llu %.1fs\n",sz,tid,(unsigned long long)ts->nproc,now()-t0);
          }
          for(int sh=0;sh<NSH;sh++){ if(ts->obn[sh]) fwrite(ts->obuf[sh],8,ts->obn[sh],ts->out[sh]);
              fclose(ts->out[sh]); }
          grand+=ts->nproc; free(ts);
        }
        fprintf(stderr,"TABLE total %llu keys zero=%llu %.1fs\n",(unsigned long long)grand,
                (unsigned long long)ZEROEV,now()-t0);
        return 0;
    }

    if(!strcmp(mode,"scan")){
        SZ=atoi(argv[3]); const char*pref=argv[4];
        for(int sh=0;sh<NSH;sh++){
            char fn[512]; snprintf(fn,sizeof fn,"%s.%d",pref,sh);
            int fd=open(fn,O_RDONLY); if(fd<0){fprintf(stderr,"missing %s\n",fn);return 2;}
            struct stat st; fstat(fd,&st); TBLN[sh]=st.st_size/8;
            TBL[sh]=(st.st_size? mmap(NULL,st.st_size,PROT_READ,MAP_SHARED,fd,0):NULL);
        }
        int fd2=open(argv[5],O_RDONLY); struct stat st2; fstat(fd2,&st2);
        BM=mmap(NULL,st2.st_size,PROT_READ,MAP_SHARED,fd2,0);
        REPORT=fopen(argv[6],"a");
        int lo=(argc>7)?atoi(argv[7]):0, hi=(argc>8)?atoi(argv[8]):512;
        static fe A1x[512],A1y[512];
        { TS*t=calloc(1,sizeof(TS)); fe ox[512],oy[512]; batch_full(t,BX,BY,0,ox,oy);
          for(int i=0;i<512;i++){fe_copy(A1x[i],ox[i]);fe_copy(A1y[i],oy[i]);} free(t); }
        u64 total=0;
        if(SZ==1){
            for(int s=lo;s<hi;s++) if(tbl_has(A1x[s][0]))
                fprintf(REPORT,"HIT 1 %d %d %llu\n",s,s,(unsigned long long)A1x[s][0]);
            total=hi-lo;
        } else {
            #pragma omp parallel reduction(+:total)
            { TS*ts=calloc(1,sizeof(TS));
              #pragma omp for schedule(dynamic,1)
              for(int s0=lo;s0<hi;s0++){
                  if(nxt(s0) > 512-2*(SZ-1)) continue;
                  rec(ts,A1x[s0],A1y[s0],nxt(s0),1,(u64)s0);
                  flush_last(ts);
              }
              flush_last(ts); total+=ts->nproc; free(ts); }
        }
        fprintf(REPORT,"DONE %s sz=%d range=[%d,%d) n=%llu zero=%llu %.1fs\n",argv[2],SZ,lo,hi,
                (unsigned long long)total,(unsigned long long)ZEROEV,now()-t0);
        fflush(REPORT);
        fprintf(stderr,"scan %s sz=%d [%d,%d): %llu cand zero=%llu %.1fs\n",argv[2],SZ,lo,hi,
                (unsigned long long)total,(unsigned long long)ZEROEV,now()-t0);
        return 0;
    }
    return 1;
}
