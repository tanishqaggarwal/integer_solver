/* agent AE -- interval Pollard-kangaroo with distinguished points, secp256k1 field.
 *
 * Solves  log_G(Q) in [0, 2^R)  for a supplied (already shifted) target Q.
 * O(1) memory in the range; memory is only the DP table.
 *
 * Usage:  aekang <inputfile> <threads> <kang_per_thread> <dpbits> <log2maxjumps> <seed> <tablebits>
 *
 * Input file (text, all big integers in hex, no 0x):
 *   line 1: Gx Gy
 *   line 2: Qx Qy          (the shifted target; log_G Q is claimed to lie in [0,2^R))
 *   line 3: R              (decimal)
 *   line 4: NJ             (decimal, power of two)
 *   next NJ lines: jd Jx Jy     (jd = jump distance, hex; J = jd*G)
 *   next 256 lines: Lx Ly       (ladder: L_i = 2^i * G)
 *
 * Output: STATUS lines, CAND lines, and a terminal DONE/FOUND line with counters.
 * Every candidate is printed only; verification is done elsewhere in bignum Python.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;
typedef struct { u64 v[4]; } fe;

#define RC 0x1000003D1ULL          /* 2^256 mod p = 2^32 + 977 */
static const u64 PP[4] = {0xFFFFFFFEFFFFFC2FULL, ~0ULL, ~0ULL, ~0ULL};

/* ------------------------------------------------------------------ field */
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
    if(c){ /* wrapped 2^256: add RC */
        u64 cc=RC;
        for(int i=0;i<4 && cc;i++){ u128 s=(u128)r->v[i]+cc; r->v[i]=(u64)s; cc=(u64)(s>>64); }
    }
    if(fe_ge_p(r)) fe_subp(r);
}
static inline void fe_sub(fe *r,const fe*a,const fe*b){
    u64 br=0;
    for(int i=0;i<4;i++){ u128 s=(u128)a->v[i]-b->v[i]-br; r->v[i]=(u64)s; br=(u64)((s>>127)&1); }
    if(br){ /* borrowed 2^256: subtract RC */
        u64 cc=RC;
        for(int i=0;i<4 && cc;i++){ u128 s=(u128)r->v[i]-cc; r->v[i]=(u64)s; cc=(u64)((s>>127)&1); }
    }
}
static inline int fe_iszero(const fe*a){ return !(a->v[0]|a->v[1]|a->v[2]|a->v[3]); }
static inline int fe_eq(const fe*a,const fe*b){
    return a->v[0]==b->v[0]&&a->v[1]==b->v[1]&&a->v[2]==b->v[2]&&a->v[3]==b->v[3];
}
static void fe_mul(fe *r,const fe*a,const fe*b){
    u64 t[8]; memset(t,0,sizeof(t));
    for(int i=0;i<4;i++){
        u64 carry=0;
        for(int j=0;j<4;j++){
            u128 s=(u128)a->v[i]*b->v[j] + t[i+j] + carry;
            t[i+j]=(u64)s; carry=(u64)(s>>64);
        }
        t[i+4]=carry;
    }
    u64 c[5],carry=0;
    for(int i=0;i<4;i++){ u128 s=(u128)t[4+i]*RC + t[i] + carry; c[i]=(u64)s; carry=(u64)(s>>64); }
    c[4]=carry;
    u128 e=(u128)c[4]*RC + c[0];
    c[0]=(u64)e; u64 cr=(u64)(e>>64);
    for(int i=1;i<4 && cr;i++){ u128 s=(u128)c[i]+cr; c[i]=(u64)s; cr=(u64)(s>>64); }
    if(cr){ u64 cc=RC; for(int i=0;i<4 && cc;i++){ u128 s=(u128)c[i]+cc; c[i]=(u64)s; cc=(u64)(s>>64);} }
    memcpy(r->v,c,32);
    if(fe_ge_p(r)) fe_subp(r);
}
static void fe_inv(fe *r,const fe*a){          /* Fermat: a^(p-2) */
    /* p-2 = 2^256 - 2^32 - 979 */
    static const u64 E[4]={0xFFFFFFFEFFFFFC2DULL,~0ULL,~0ULL,~0ULL};
    fe x=*a, acc; int started=0;
    memset(&acc,0,sizeof(acc)); acc.v[0]=1;
    for(int i=255;i>=0;i--){
        if(started) fe_mul(&acc,&acc,&acc);
        if((E[i>>6]>>(i&63))&1){ if(!started){ acc=x; started=1; } else fe_mul(&acc,&acc,&x); }
    }
    *r=acc;
}

/* ------------------------------------------------------------------ points */
typedef struct { fe x,y; int inf; } pt;

static void pt_add(pt *r,const pt*P,const pt*Q){   /* generic affine add, own inversion */
    if(P->inf){ *r=*Q; return; } if(Q->inf){ *r=*P; return; }
    fe dx,dy,l,t1,t2;
    fe_sub(&dx,&Q->x,&P->x);
    if(fe_iszero(&dx)){
        fe s; fe_add(&s,&P->y,&Q->y);
        if(fe_iszero(&s)){ r->inf=1; return; }
        fe_mul(&t1,&P->x,&P->x); fe_add(&t2,&t1,&t1); fe_add(&t2,&t2,&t1); /* 3x^2 */
        fe_add(&dy,&P->y,&P->y); fe_inv(&dx,&dy); fe_mul(&l,&t2,&dx);
    } else {
        fe_sub(&dy,&Q->y,&P->y); fe_inv(&t1,&dx); fe_mul(&l,&dy,&t1);
    }
    fe_mul(&t1,&l,&l); fe_sub(&t1,&t1,&P->x); fe_sub(&t1,&t1,&Q->x);
    fe_sub(&t2,&P->x,&t1); fe_mul(&t2,&l,&t2); fe_sub(&t2,&t2,&P->y);
    r->x=t1; r->y=t2; r->inf=0;
}

/* ------------------------------------------------------------------ globals */
static pt Gp, Qp, Jp[64], Lad[256];
static u128 Jd[64];
static int NJ, Rbits;
static int NTH, KPT, DPBITS, TBITS;
static u64 MAXJ;
static u64 DPMASK;

typedef struct { u64 x0,x1; u128 s; uint8_t used,type; } dpent;
static dpent *tab; static u64 TSZ, TMASK;
static pthread_mutex_t tlock = PTHREAD_MUTEX_INITIALIZER;
static volatile u64 g_jumps=0, g_dps=0, g_used=0, g_cands=0, g_dxzero=0;
static volatile int g_stop=0;
static pthread_mutex_t clock_ = PTHREAD_MUTEX_INITIALIZER;

static void print_u128(char*buf,u128 v){
    char tmp[48]; int n=0;
    if(v==0){ strcpy(buf,"0"); return; }
    while(v){ tmp[n++]='0'+(int)(v%10); v/=10; }
    for(int i=0;i<n;i++) buf[i]=tmp[n-1-i];
    buf[n]=0;
}

/* xorshift rng */
static inline u64 xrand(u64*s){ u64 x=*s; x^=x<<13; x^=x>>7; x^=x<<17; *s=x; return x; }

static void scalar_mul_ladder(pt*r,u128 k){      /* k*G via the supplied 2^i ladder */
    r->inf=1;
    for(int i=0;i<128;i++) if((k>>i)&1) pt_add(r,r,&Lad[i]);
}

/* insert / probe DP table.  returns 1 if a tame-wild collision was reported */
static void dp_report(const fe*x,u128 s,int type){
    u64 k0=x->v[0], k1=x->v[1];
    u64 h=(k0*0x9E3779B97F4A7C15ULL)^(k1*0xC2B2AE3D27D4EB4FULL);
    h&=TMASK;
    pthread_mutex_lock(&tlock);
    for(u64 i=0;i<TSZ;i++){
        u64 idx=(h+i)&TMASK;
        if(!tab[idx].used){
            tab[idx].used=1; tab[idx].x0=k0; tab[idx].x1=k1; tab[idx].s=s; tab[idx].type=(uint8_t)type;
            g_used++;
            break;
        }
        if(tab[idx].x0==k0 && tab[idx].x1==k1){
            if(tab[idx].type!=type){
                char a[48],b[48];
                u128 st = (type==0)? s : tab[idx].s;
                u128 sw = (type==0)? tab[idx].s : s;
                print_u128(a,st); print_u128(b,sw);
                printf("CAND st=%s sw=%s\n",a,b); fflush(stdout);
                g_cands++; g_stop=1;   /* any tame-wild DP collision is an exact relation */
            }
            break;
        }
    }
    if(g_used*10 > TSZ*6){ fprintf(stderr,"DPTABLE_FULL\n"); g_stop=1; }
    pthread_mutex_unlock(&tlock);
}

typedef struct { int id; u64 seed; } targ;

static void* worker(void*vp){
    targ*ta=(targ*)vp;
    int K=KPT;
    u64 rs = ta->seed ? ta->seed : 0x243F6A8885A308D3ULL;
    fe *X=malloc(sizeof(fe)*K), *Y=malloc(sizeof(fe)*K);
    u128 *S=malloc(sizeof(u128)*K);
    uint8_t *TY=malloc(K);
    u64 *since=calloc(K,sizeof(u64));
    fe *dxs=malloc(sizeof(fe)*K), *pref=malloc(sizeof(fe)*K);
    u128 L = (Rbits>=128)? (u128)0 : ((u128)1<<Rbits);
    u128 Wwin = L>>3; if(Wwin==0) Wwin=1;

    for(int i=0;i<K;i++){
        int tame = (i&1)==0;
        u128 off = ((u128)xrand(&rs)<<64 | xrand(&rs));
        pt P;
        if(tame){ off %= L; scalar_mul_ladder(&P,off); }
        else { off %= Wwin; pt A; scalar_mul_ladder(&A,off); pt_add(&P,&Qp,&A); }
        X[i]=P.x; Y[i]=P.y; S[i]=off; TY[i]=tame?0:1;
    }

    u64 local=0;
    while(!g_stop){
        /* one batched step over all K kangaroos */
        for(int i=0;i<K;i++){
            int ji = (int)(X[i].v[0] & (NJ-1));
            fe_sub(&dxs[i], &Jp[ji].x, &X[i]);
        }
        /* Montgomery batch inversion */
        fe run; memset(&run,0,sizeof(run)); run.v[0]=1;
        int zero_seen=0;
        for(int i=0;i<K;i++){
            pref[i]=run;
            if(fe_iszero(&dxs[i])){ zero_seen=1; }
            else fe_mul(&run,&run,&dxs[i]);
        }
        fe invrun; fe_inv(&invrun,&run);
        for(int i=K-1;i>=0;i--){
            if(fe_iszero(&dxs[i])) continue;
            fe t; fe_mul(&t,&invrun,&pref[i]);       /* = 1/dxs[i] */
            fe_mul(&invrun,&invrun,&dxs[i]);
            dxs[i]=t;
        }
        for(int i=0;i<K;i++){
            int ji=(int)(X[i].v[0] & (NJ-1));
            if(fe_iszero(&dxs[i])){                  /* degenerate: re-randomise */
                __sync_fetch_and_add(&g_dxzero,1);
                u128 off=((u128)xrand(&rs)<<64|xrand(&rs));
                pt P;
                if(TY[i]==0){ off%=L; scalar_mul_ladder(&P,off); }
                else { off%=Wwin; pt A; scalar_mul_ladder(&A,off); pt_add(&P,&Qp,&A); }
                X[i]=P.x; Y[i]=P.y; S[i]=off; since[i]=0; continue;
            }
            fe dy,l,nx,ny,t;
            fe_sub(&dy,&Jp[ji].y,&Y[i]);
            fe_mul(&l,&dy,&dxs[i]);
            fe_mul(&nx,&l,&l); fe_sub(&nx,&nx,&X[i]); fe_sub(&nx,&nx,&Jp[ji].x);
            fe_sub(&t,&X[i],&nx); fe_mul(&ny,&l,&t); fe_sub(&ny,&ny,&Y[i]);
            X[i]=nx; Y[i]=ny; S[i]+=Jd[ji]; since[i]++;
            if((X[i].v[0]&DPMASK)==0){
                dp_report(&X[i],S[i],TY[i]);
                __sync_fetch_and_add(&g_dps,1);
                since[i]=0;
            } else if(since[i] > (u64)40<<DPBITS){    /* trapped: re-randomise */
                u128 off=((u128)xrand(&rs)<<64|xrand(&rs));
                pt P;
                if(TY[i]==0){ off%=L; scalar_mul_ladder(&P,off); }
                else { off%=Wwin; pt A; scalar_mul_ladder(&A,off); pt_add(&P,&Qp,&A); }
                X[i]=P.x; Y[i]=P.y; S[i]=off; since[i]=0;
            }
        }
        (void)zero_seen;
        local+=K;
        if(local>=(u64)K*64){
            u64 tot=__sync_add_and_fetch(&g_jumps,local); local=0;
            if(tot>=MAXJ) g_stop=1;
        }
    }
    __sync_add_and_fetch(&g_jumps,local);
    free(X);free(Y);free(S);free(TY);free(since);free(dxs);free(pref);
    return NULL;
}

static int rdfe(FILE*f,fe*a){
    char buf[80];
    if(fscanf(f,"%79s",buf)!=1) return 0;
    memset(a,0,sizeof(fe));
    int n=strlen(buf);
    for(int i=0;i<n;i++){
        int c=buf[i], d;
        if(c>='0'&&c<='9') d=c-'0'; else if(c>='a'&&c<='f') d=c-'a'+10;
        else if(c>='A'&&c<='F') d=c-'A'+10; else return 0;
        /* a = a*16 + d */
        u64 carry=d;
        for(int k=0;k<4;k++){ u64 nv=(a->v[k]<<4)|carry; carry=a->v[k]>>60; a->v[k]=nv; }
    }
    return 1;
}

int main(int argc,char**argv){
    if(argc<8){ fprintf(stderr,"usage: %s in threads kpt dpbits log2maxjumps seed tablebits\n",argv[0]); return 2; }
    FILE*f=fopen(argv[1],"r"); if(!f){ perror("open"); return 2; }
    NTH=atoi(argv[2]); KPT=atoi(argv[3]); DPBITS=atoi(argv[4]);
    int lgmax=atoi(argv[5]); u64 seed=strtoull(argv[6],NULL,10); TBITS=atoi(argv[7]);
    MAXJ = (lgmax>=63)? (u64)1<<62 : ((u64)1<<lgmax);
    DPMASK = (DPBITS>=64)?~0ULL:(((u64)1<<DPBITS)-1);

    if(!rdfe(f,&Gp.x)||!rdfe(f,&Gp.y)) { fprintf(stderr,"bad G\n"); return 2; } Gp.inf=0;
    if(!rdfe(f,&Qp.x)||!rdfe(f,&Qp.y)) { fprintf(stderr,"bad Q\n"); return 2; } Qp.inf=0;
    if(fscanf(f,"%d",&Rbits)!=1){ fprintf(stderr,"bad R\n"); return 2; }
    if(fscanf(f,"%d",&NJ)!=1||NJ>64||(NJ&(NJ-1))){ fprintf(stderr,"bad NJ\n"); return 2; }
    for(int i=0;i<NJ;i++){
        fe t; if(!rdfe(f,&t)) { fprintf(stderr,"bad jd\n"); return 2; }
        Jd[i]=((u128)t.v[1]<<64)|t.v[0];
        if(!rdfe(f,&Jp[i].x)||!rdfe(f,&Jp[i].y)){ fprintf(stderr,"bad J\n"); return 2; } Jp[i].inf=0;
    }
    for(int i=0;i<256;i++){
        if(!rdfe(f,&Lad[i].x)||!rdfe(f,&Lad[i].y)){ fprintf(stderr,"bad ladder %d\n",i); return 2; }
        Lad[i].inf=0;
    }
    fclose(f);

    /* self-test: ladder consistency (L_{i+1} = 2 L_i) on a few indices */
    for(int i=0;i<8;i++){
        pt d; pt_add(&d,&Lad[i],&Lad[i]);
        if(!fe_eq(&d.x,&Lad[i+1].x)||!fe_eq(&d.y,&Lad[i+1].y)){
            fprintf(stderr,"SELFTEST_FAIL ladder %d\n",i); return 3; }
    }
    /* self-test: J_i == jd_i * G */
    for(int i=0;i<NJ;i++){
        pt c; scalar_mul_ladder(&c,Jd[i]);
        if(!fe_eq(&c.x,&Jp[i].x)||!fe_eq(&c.y,&Jp[i].y)){
            fprintf(stderr,"SELFTEST_FAIL jump %d\n",i); return 3; }
    }
    fprintf(stderr,"SELFTEST_OK\n");

    TSZ=(u64)1<<TBITS; TMASK=TSZ-1;
    tab=calloc(TSZ,sizeof(dpent));
    if(!tab){ fprintf(stderr,"table alloc failed\n"); return 2; }

    pthread_t*th=malloc(sizeof(pthread_t)*NTH);
    targ*ta=malloc(sizeof(targ)*NTH);
    struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
    for(int i=0;i<NTH;i++){ ta[i].id=i; ta[i].seed=seed+0x9E3779B97F4A7C15ULL*(i+1); pthread_create(&th[i],NULL,worker,&ta[i]); }
    /* monitor */
    int tick=0;
    while(!g_stop){
        struct timespec ts={1,0}; nanosleep(&ts,NULL);
        clock_gettime(CLOCK_MONOTONIC,&t1);
        double el=(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
        if(++tick%10==0){
            fprintf(stderr,"STATUS jumps=%llu dps=%llu cands=%llu t=%.1f rate=%.3fM/s\n",
                (unsigned long long)g_jumps,(unsigned long long)g_dps,
                (unsigned long long)g_cands,el,g_jumps/el/1e6);
        }
        if(g_jumps>=MAXJ) g_stop=1;
    }
    for(int i=0;i<NTH;i++) pthread_join(th[i],NULL);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    double el=(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
    printf("DONE jumps=%llu dps=%llu dpexp=%.1f cands=%llu dxzero=%llu tabused=%llu secs=%.1f rate=%.3fM/s\n",
        (unsigned long long)g_jumps,(unsigned long long)g_dps,
        (double)g_jumps/(double)((u64)1<<DPBITS),
        (unsigned long long)g_cands,(unsigned long long)g_dxzero,
        (unsigned long long)g_used,el,g_jumps/el/1e6);
    return 0;
}
