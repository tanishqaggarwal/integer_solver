/* agent AE -- small-quotient sweep.
 *
 * Finds all (a,b,e) with 1 <= a,b <= M and  a*G == +/- lam^e * (b*T),
 * i.e.  k0 == +/- a * lam^-e * b^-1  (mod N).
 *
 * Two incremental walks (P += G and P += T) run in W lanes with one batched
 * inversion per lane-step; the G-side x-coordinates are stored and sorted, the
 * T-side x-coordinates (and their beta, beta^2 multiples, which are the
 * x-coordinates of the lam and lam^2 images) are looked up.
 *
 * usage: aequot <input> <out>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;
typedef struct { u64 v[4]; } fe;
#define RC 0x1000003D1ULL
static const u64 PP[4] = {0xFFFFFFFEFFFFFC2FULL, ~0ULL, ~0ULL, ~0ULL};
static inline int fe_ge_p(const fe *a){
    for(int i=3;i>=0;i--){ if(a->v[i]>PP[i]) return 1; if(a->v[i]<PP[i]) return 0; }
    return 1;
}
static inline void fe_subp(fe *a){
    u128 br=0; u64 r[4];
    for(int i=0;i<4;i++){ u128 s=(u128)a->v[i]-PP[i]-br; r[i]=(u64)s; br=(s>>127)&1; }
    if(!br) memcpy(a->v,r,32);
}
static inline void fe_add(fe *r,const fe*a,const fe*b){
    u64 c=0;
    for(int i=0;i<4;i++){ u128 s=(u128)a->v[i]+b->v[i]+c; r->v[i]=(u64)s; c=(u64)(s>>64); }
    if(c){ u64 cc=RC; for(int i=0;i<4&&cc;i++){ u128 s=(u128)r->v[i]+cc; r->v[i]=(u64)s; cc=(u64)(s>>64);} }
    if(fe_ge_p(r)) fe_subp(r);
}
static inline void fe_sub(fe *r,const fe*a,const fe*b){
    u64 br=0;
    for(int i=0;i<4;i++){ u128 s=(u128)a->v[i]-b->v[i]-br; r->v[i]=(u64)s; br=(u64)((s>>127)&1); }
    if(br){ u64 cc=RC; for(int i=0;i<4&&cc;i++){ u128 s=(u128)r->v[i]-cc; r->v[i]=(u64)s; cc=(u64)((s>>127)&1);} }
}
static inline int fe_iszero(const fe*a){ return !(a->v[0]|a->v[1]|a->v[2]|a->v[3]); }
static void fe_mul(fe *r,const fe*a,const fe*b){
    u64 t[8]; memset(t,0,sizeof(t));
    for(int i=0;i<4;i++){
        u64 carry=0;
        for(int j=0;j<4;j++){ u128 s=(u128)a->v[i]*b->v[j]+t[i+j]+carry; t[i+j]=(u64)s; carry=(u64)(s>>64); }
        t[i+4]=carry;
    }
    u64 c[5],carry=0;
    for(int i=0;i<4;i++){ u128 s=(u128)t[4+i]*RC+t[i]+carry; c[i]=(u64)s; carry=(u64)(s>>64); }
    c[4]=carry;
    u128 e=(u128)c[4]*RC+c[0]; c[0]=(u64)e; u64 cr=(u64)(e>>64);
    for(int i=1;i<4&&cr;i++){ u128 s=(u128)c[i]+cr; c[i]=(u64)s; cr=(u64)(s>>64); }
    if(cr){ u64 cc=RC; for(int i=0;i<4&&cc;i++){ u128 s=(u128)c[i]+cc; c[i]=(u64)s; cc=(u64)(s>>64);} }
    memcpy(r->v,c,32);
    if(fe_ge_p(r)) fe_subp(r);
}
static void fe_inv(fe *r,const fe*a){
    static const u64 E[4]={0xFFFFFFFEFFFFFC2DULL,~0ULL,~0ULL,~0ULL};
    fe x=*a,acc; int started=0; memset(&acc,0,sizeof(acc)); acc.v[0]=1;
    for(int i=255;i>=0;i--){
        if(started) fe_mul(&acc,&acc,&acc);
        if((E[i>>6]>>(i&63))&1){ if(!started){acc=x;started=1;} else fe_mul(&acc,&acc,&x); }
    }
    *r=acc;
}
static int rdfe(FILE*f,fe*a){
    char buf[80]; if(fscanf(f,"%79s",buf)!=1) return 0;
    memset(a,0,sizeof(fe)); int n=strlen(buf);
    for(int i=0;i<n;i++){
        int c=buf[i],d;
        if(c>='0'&&c<='9')d=c-'0'; else if(c>='a'&&c<='f')d=c-'a'+10;
        else if(c>='A'&&c<='F')d=c-'A'+10; else return 0;
        u64 carry=d;
        for(int k=0;k<4;k++){ u64 nv=(a->v[k]<<4)|carry; carry=a->v[k]>>60; a->v[k]=nv; }
    }
    return 1;
}

typedef struct { u64 key; u64 val; } ent;
static int cmpent(const void*A,const void*B){
    u64 x=((const ent*)A)->key, y=((const ent*)B)->key;
    return (x<y)?-1:((x>y)?1:0);
}

int main(int argc,char**argv){
    if(argc<3){ fprintf(stderr,"usage: %s in out\n",argv[0]); return 2; }
    FILE*f=fopen(argv[1],"r"); if(!f){perror("in");return 2;}
    fe Gx,Gy,Tx,Ty,BETA;
    int LOGM,W;
    if(!rdfe(f,&Gx)||!rdfe(f,&Gy)||!rdfe(f,&Tx)||!rdfe(f,&Ty)||!rdfe(f,&BETA)){fprintf(stderr,"bad hdr\n");return 2;}
    if(fscanf(f,"%d %d",&LOGM,&W)!=2){fprintf(stderr,"bad LOGM/W\n");return 2;}
    u64 M=(u64)1<<LOGM;
    if(M % (u64)W){ fprintf(stderr,"W must divide M\n"); return 2; }
    u64 chunk=M/W;
    fe *AX=malloc(sizeof(fe)*W),*AY=malloc(sizeof(fe)*W);
    fe *BX=malloc(sizeof(fe)*W),*BY=malloc(sizeof(fe)*W);
    for(int j=0;j<W;j++) if(!rdfe(f,&AX[j])||!rdfe(f,&AY[j])){fprintf(stderr,"bad A%d\n",j);return 2;}
    for(int j=0;j<W;j++) if(!rdfe(f,&BX[j])||!rdfe(f,&BY[j])){fprintf(stderr,"bad B%d\n",j);return 2;}
    fclose(f);

    ent *arr=malloc(sizeof(ent)*M);
    if(!arr){ fprintf(stderr,"alloc %llu MB failed\n",(unsigned long long)(sizeof(ent)*M>>20)); return 2; }
    fe *dxs=malloc(sizeof(fe)*W),*pref=malloc(sizeof(fe)*W);
    FILE*o=fopen(argv[2],"w");
    struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);

    /* ---- phase 1: a*G for a = j*chunk+1 .. (j+1)*chunk ---- */
    u64 written=0, deg1=0;
    for(u64 s=0;s<chunk;s++){
        for(int j=0;j<W;j++){ arr[written].key=AX[j].v[0]; arr[written].val=(u64)j*chunk+1+s; written++; }
        if(s+1==chunk) break;
        for(int j=0;j<W;j++) fe_sub(&dxs[j],&Gx,&AX[j]);
        fe run; memset(&run,0,sizeof(run)); run.v[0]=1;
        for(int j=0;j<W;j++){ pref[j]=run; if(!fe_iszero(&dxs[j])) fe_mul(&run,&run,&dxs[j]); }
        fe ir; fe_inv(&ir,&run);
        for(int j=W-1;j>=0;j--){ if(fe_iszero(&dxs[j])) continue; fe t; fe_mul(&t,&ir,&pref[j]); fe_mul(&ir,&ir,&dxs[j]); dxs[j]=t; }
        for(int j=0;j<W;j++){
            if(fe_iszero(&dxs[j])){ deg1++; continue; }
            fe dy,l,nx,ny,t;
            fe_sub(&dy,&Gy,&AY[j]); fe_mul(&l,&dy,&dxs[j]);
            fe_mul(&nx,&l,&l); fe_sub(&nx,&nx,&AX[j]); fe_sub(&nx,&nx,&Gx);
            fe_sub(&t,&AX[j],&nx); fe_mul(&ny,&l,&t); fe_sub(&ny,&ny,&AY[j]);
            AX[j]=nx; AY[j]=ny;
        }
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    fprintf(stderr,"PHASE1 written=%llu expect=%llu deg=%llu %.1fs\n",
        (unsigned long long)written,(unsigned long long)M,(unsigned long long)deg1,
        (t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec));
    if(written!=M){ fprintf(stderr,"PHASE1_COUNT_MISMATCH\n"); return 3; }

    qsort(arr,M,sizeof(ent),cmpent);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    fprintf(stderr,"SORTED %.1fs\n",(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec));

    /* ---- phase 2: b*T, probe x, beta*x, beta^2*x ---- */
    fe B2; fe_mul(&B2,&BETA,&BETA);
    u64 probed=0, hits=0, deg2=0;
    for(u64 s=0;s<chunk;s++){
        for(int j=0;j<W;j++){
            u64 bb=(u64)j*chunk+1+s;
            fe q=BX[j];
            for(int e=0;e<3;e++){
                if(e==1) fe_mul(&q,&BX[j],&BETA);
                if(e==2) fe_mul(&q,&BX[j],&B2);
                u64 key=q.v[0];
                /* binary search for key */
                u64 lo=0,hi=M;
                while(lo<hi){ u64 mid=(lo+hi)>>1; if(arr[mid].key<key) lo=mid+1; else hi=mid; }
                while(lo<M && arr[lo].key==key){
                    fprintf(o,"MATCH a=%llu b=%llu e=%d\n",(unsigned long long)arr[lo].val,
                            (unsigned long long)bb,e);
                    hits++; lo++;
                }
                probed++;
            }
        }
        if(s+1==chunk) break;
        for(int j=0;j<W;j++) fe_sub(&dxs[j],&Tx,&BX[j]);
        fe run; memset(&run,0,sizeof(run)); run.v[0]=1;
        for(int j=0;j<W;j++){ pref[j]=run; if(!fe_iszero(&dxs[j])) fe_mul(&run,&run,&dxs[j]); }
        fe ir; fe_inv(&ir,&run);
        for(int j=W-1;j>=0;j--){ if(fe_iszero(&dxs[j])) continue; fe t; fe_mul(&t,&ir,&pref[j]); fe_mul(&ir,&ir,&dxs[j]); dxs[j]=t; }
        for(int j=0;j<W;j++){
            if(fe_iszero(&dxs[j])){ deg2++; continue; }
            fe dy,l,nx,ny,t;
            fe_sub(&dy,&Ty,&BY[j]); fe_mul(&l,&dy,&dxs[j]);
            fe_mul(&nx,&l,&l); fe_sub(&nx,&nx,&BX[j]); fe_sub(&nx,&nx,&Tx);
            fe_sub(&t,&BX[j],&nx); fe_mul(&ny,&l,&t); fe_sub(&ny,&ny,&BY[j]);
            BX[j]=nx; BY[j]=ny;
        }
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    double el=(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
    fprintf(o,"DONE M=%llu probes=%llu expect=%llu matches=%llu deg1=%llu deg2=%llu secs=%.1f\n",
        (unsigned long long)M,(unsigned long long)probed,(unsigned long long)(3*M),
        (unsigned long long)hits,(unsigned long long)deg1,(unsigned long long)deg2,el);
    fprintf(stderr,"DONE probes=%llu expect=%llu matches=%llu %.1fs\n",
        (unsigned long long)probed,(unsigned long long)(3*M),(unsigned long long)hits,el);
    fclose(o);
    return (probed==3*M)?0:3;
}
