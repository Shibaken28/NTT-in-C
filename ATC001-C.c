# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define S 262144
# define LOG_S 18
# define MOD 998244353

# define DIVIDE_LIMIT 23

long primitive_root=3;
long primitive_pow[DIVIDE_LIMIT+1];//[i] := i^2乗根
long primitive_pow_inv[DIVIDE_LIMIT+1];// [i] := i^2乗根の逆数

long modPow(long a,long b,long m){
    if(b==0)return 1;
    if(b==1)return a%m;
    if(b%2==1)return (a*modPow(a,b-1,m))%m;
    long c=modPow(a,b/2,m);
    return (c*c)%m;
}
long modInv(long a,long m){
	return modPow(a,m-2,m);
}



void calmod(){
    primitive_pow[DIVIDE_LIMIT]=modPow(primitive_root,119,MOD);//3^119 = 2^DIVIDE_LIMIT乗根
    //15311432
    primitive_pow_inv[DIVIDE_LIMIT]=modInv(primitive_pow[DIVIDE_LIMIT],MOD);//それの逆数
    //469870224
    for(int i=DIVIDE_LIMIT-1;i>=0;i--){
        primitive_pow[i] = (primitive_pow[i+1]*primitive_pow[i+1])%MOD;
        primitive_pow_inv[i] = (primitive_pow_inv[i+1]*primitive_pow_inv[i+1])%MOD;
        //printf("pp[%2d] = %ld\n",i,primitive_pow[i+1]);
        //printf("pi[%2d] = %ld\n",i,primitive_pow_inv[i+1]);
    }
}


void fft(long f[],long N,long k,long r[],int inv){
    if(N==1){//1乗根...x=1
        //printf("%ld\n",f[0]);
        r[0]=f[0];
        //printf("%ld\n",r[0]);
        return;
    }
    
    long *f1,*f2;
    long *r1,*r2;
    f1 = (long*)malloc(sizeof(long)*(N/2));
    f2 = (long*)malloc(sizeof(long)*(N/2));
    r1 = (long*)malloc(sizeof(long)*(N));
    r2 = (long*)malloc(sizeof(long)*(N));
    // f(x) = f1(x^2) + x f2(x^2)
    // それぞれの結果をr1,r2
    for(int i=0;i<N;i++){//係数をコピー
        if(i%2==0){
            f1[i/2]=f[i];
        }else{
            f2[i/2]=f[i];
        }
    }
    //0,2,4,6,...2(N-1) で2周しているので前半と後半の値は同じ
    //半分だけ求めてもらう
    fft(f1,N/2,k-1,r1,inv);
    fft(f2,N/2,k-1,r2,inv);

    //半分なので2倍する
    for(int i=0;i<N/2;i++){
        r1[i+N/2]=r1[i];
        r2[i+N/2]=r2[i];
    }

    //printf("N = %d\n",N);
    for(long i=0;i<N;i++){
        //r[i]=r1[i] + theta*i r2[i];
        r[i] = r1[i];
        long x;
        if(inv==1)x = modPow(primitive_pow[k],i,MOD);
        if(inv==-1)x = modPow(primitive_pow_inv[k],i,MOD);
        r[i] += (x*r2[i]) %MOD;
        r[i]%=MOD;
    }
    free(f1);
    free(f2);
    free(r1);
    free(r2);

}


long ff[S];
long gg[S];
void def(long f[],long g[],long r[]){
    long N=S;
    fft(f,N,LOG_S,ff,1);
    fft(g,N,LOG_S,gg,1);
    
    for(int i=0;i<N;i++){
        ff[i]=(gg[i]*ff[i])%MOD;
    }
    fft(ff,N,LOG_S,r,-1);
}

long A[S];
long B[S];
long C[S];

int main(){    
    int N;
    
    scanf("%d",&N);
    for(int i=0;i<S;i++){
        A[i]=0;B[i]=0;
    }
    for(int i=1;i<=N;i++){
        scanf("%ld",&(A[i]));
        scanf("%ld",&(B[i]));
    }
    calmod();
    def(A,B,C);
    for(int i=1;i<=2*N;i++){
        printf("%ld\n",(C[i]*modInv(S,MOD))%MOD);
    }
    
    return 0;
}
