/* Agent X: baby-step giant-step for a SMALL scalar, C engine.  Covers k <= 2^52 and N-k <= 2^52.
   modes:  baby <data> <out>      -> BABY keys, x(j*G) for j=1..BABY
           giant <data> <sortedbaby> <sign> -> scan T -/+ i*BIG*G, report matches           */
#include "xfield.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#define LB 26
#define BABY (1ULL<<LB)
#define GIANT (1ULL<<LB)
#define W 512
static fe LX[256],LY[256],BX,BY;
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void hexto(u64*r,const char*s){r[0]=r[1]=r[2]=r[3]=0;
    for(const char*q=s;*q;q++){u64 c=0;for(int i=0;i<4;i++){u128 v=(u128)r[i]*10+c;r[i]=(u64)v;c=(u64)(v>>64);}
        c=*q-'0';for(int i=0;i<4&&c;i++){u128 v=(u128)r[i]+c;r[i]=(u64)v;c=(u64)(v>>64);} } }
/* generic affine add of W independent pairs: A_b += B_b  (batch inverted) */
static fe dd[W],pr[W];
static void badd(fe*ax,fe*ay,const fe*bx,const fe*by,int n,int*degen){
    for(int i=0;i<n;i++){ fe_sub(dd[i],ax[i],bx[i]); if(fe_iszero(dd[i])){*degen=1;dd[i][0]=1;dd[i][1]=dd[i][2]=dd[i][3]=0;} }
    fe acc={1,0,0,0};
    for(int i=0;i<n;i++){fe_copy(pr[i],acc);fe_mul(acc,acc,dd[i]);}
    fe ainv; fe_inv(ainv,acc);
    for(int i=n-1;i>=0;i--){
        fe iv;fe_mul(iv,ainv,pr[i]);fe_mul(ainv,ainv,dd[i]);
        fe lam,num,x3,y3,t;
        fe_sub(num,ay[i],by[i]); fe_mul(lam,num,iv);
        fe_mul(x3,lam,lam); fe_sub(x3,x3,ax[i]); fe_sub(x3,x3,bx[i]);
        fe_sub(t,ax[i],x3); fe_mul(y3,lam,t); fe_sub(y3,y3,ay[i]);
        fe_copy(ax[i],x3); fe_copy(ay[i],y3);
    }
}
/* scalar multiply by small u64 using the precomputed ladder LX/LY (positive) */
static void smul(fe rx,fe ry,u64 k,int*ok){
    int first=1; fe ax[1],ay[1],bx[1],by[1]; int dg=0;
    for(int i=0;i<64;i++) if((k>>i)&1){
        if(first){fe_copy(rx,LX[i]);fe_copy(ry,LY[i]);first=0;}
        else{fe_copy(ax[0],rx);fe_copy(ay[0],ry);fe_copy(bx[0],LX[i]);fe_copy(by[0],LY[i]);
             badd(ax,ay,bx,by,1,&dg); fe_copy(rx,ax[0]);fe_copy(ry,ay[0]);}
    }
    *ok=!first && !dg;
}
int main(int argc,char**argv){
    FILE*f=fopen(argv[2],"r"); char a[200],b[200];
    if(fscanf(f,"%199s %199s",a,b)!=2)return 1; hexto(BX,a);hexto(BY,b);
    for(int i=0;i<256;i++){ if(fscanf(f,"%199s %199s",a,b)!=2)return 1; hexto(LX[i],a);hexto(LY[i],b);} fclose(f);
    double t0=now();
    if(!strcmp(argv[1],"baby")){
        /* W chains: chain c holds (c+1 + m*W)*G ; step by W*G */
        static fe cx[W],cy[W],sx[W],sy[W]; int ok;
        fe SX,SY; smul(SX,SY,(u64)W,&ok); if(!ok){fprintf(stderr,"smul fail\n");return 1;}
        for(int c=0;c<W;c++){ smul(cx[c],cy[c],(u64)(c+1),&ok); if(!ok){fprintf(stderr,"init fail %d\n",c);return 1;}
                              fe_copy(sx[c],SX); fe_copy(sy[c],SY); }
        FILE*o=fopen(argv[3],"wb");
        u64 *ob=malloc(sizeof(u64)*W); u64 n=0; int dg=0;
        u64 steps=BABY/W;
        for(u64 m=0;m<steps;m++){
            for(int c=0;c<W;c++) ob[c]=cx[c][0];
            fwrite(ob,8,W,o); n+=W;
            if(m+1<steps) badd(cx,cy,sx,sy,W,&dg);
            if((m&0xFFFF)==0) fprintf(stderr,"  baby %llu %.1fs\n",(unsigned long long)n,now()-t0);
        }
        fclose(o); fprintf(stderr,"baby %llu keys degen=%d %.1fs\n",(unsigned long long)n,dg,now()-t0);
        return 0;
    }
    if(!strcmp(argv[1],"giant")){
        int fd=open(argv[3],O_RDONLY); struct stat st; fstat(fd,&st); size_t nb=st.st_size/8;
        const u64*bt=mmap(NULL,st.st_size,PROT_READ,MAP_SHARED,fd,0);
        int sign=atoi(argv[4]);
        /* base = T (sign +1) or -T (sign -1);  Q_i = base - i*BIG*G , BIG = BABY */
        fe BIGX,BIGY; int ok; smul(BIGX,BIGY,BABY,&ok); if(!ok)return 1;
        fe_neg(BIGY,BIGY);                                  /* subtract */
        static fe cx[W],cy[W],sx[W],sy[W];
        /* chain c starts at base - c*BIG , steps by -W*BIG */
        fe bx0,by0; fe_copy(bx0,BX); fe_copy(by0,BY); if(sign<0) fe_neg(by0,by0);
        fe stepx[1],stepy[1],accx[1],accy[1]; int dg=0;
        fe_copy(accx[0],bx0); fe_copy(accy[0],by0);
        for(int c=0;c<W;c++){ fe_copy(cx[c],accx[0]); fe_copy(cy[c],accy[0]);
            fe_copy(stepx[0],BIGX); fe_copy(stepy[0],BIGY);
            badd(accx,accy,stepx,stepy,1,&dg); }
        /* per-chain step = -W*BIG : accumulate by repeatedly adding -BIG W times to a point */
        fe WSX,WSY; { fe px[1],py[1],qx[1],qy[1]; fe_copy(px[0],BIGX);fe_copy(py[0],BIGY);
            for(int t=1;t<W;t++){fe_copy(qx[0],BIGX);fe_copy(qy[0],BIGY);badd(px,py,qx,qy,1,&dg);}
            fe_copy(WSX,px[0]); fe_copy(WSY,py[0]); }
        for(int c=0;c<W;c++){fe_copy(sx[c],WSX);fe_copy(sy[c],WSY);}
        u64 steps=GIANT/W, checked=0, cand=0;
        for(u64 m=0;m<steps;m++){
            for(int c=0;c<W;c++){
                u64 key=cx[c][0];
                /* binary search */
                size_t lo=0,hi=nb; while(lo<hi){size_t mid=(lo+hi)>>1; if(bt[mid]<key)lo=mid+1; else hi=mid;}
                if(lo<nb&&bt[lo]==key){ cand++; printf("CAND sign=%d i=%llu key=%llu\n",sign,
                        (unsigned long long)(m*W+c),(unsigned long long)key); }
            }
            checked+=W;
            if(m+1<steps) badd(cx,cy,sx,sy,W,&dg);
            if((m&0xFFFF)==0) fprintf(stderr,"  giant %llu %.1fs\n",(unsigned long long)checked,now()-t0);
        }
        fprintf(stderr,"giant sign=%d checked=%llu candidates=%llu degen=%d %.1fs\n",
                sign,(unsigned long long)checked,(unsigned long long)cand,dg,now()-t0);
        printf("DONE sign=%d checked=%llu cand=%llu (covers |k| <= 2^%d)\n",sign,
               (unsigned long long)checked,(unsigned long long)cand,2*LB);
        return 0;
    }
    if(!strcmp(argv[1],"smultest")){
        for(int t=0;t<5;t++){ u64 k=(u64)1<<(10*t+3); k+=12345*t+7; int ok; fe rx,ry; smul(rx,ry,k,&ok);
            printf("smul %llu ok=%d x=%llu %llu %llu %llu\n",(unsigned long long)k,ok,
              (unsigned long long)rx[0],(unsigned long long)rx[1],(unsigned long long)rx[2],(unsigned long long)rx[3]); }
        return 0; }
    return 1;
}
