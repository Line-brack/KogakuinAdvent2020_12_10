//参考
//http://sys.ci.ritsumei.ac.jp/probability/2011/10-6.pdf
//https://nishiru3.hatenablog.com/entry/2014/10/28/085130
//https://qiita.com/piyo7/items/36b19c6ea68baa2a69cd
#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;
const double pi=acos(-1);
//合同式法による一様な疑似乱数生成
class CongruenceMethod{
    private:
    int a,c,m;
    int r;
    double x;
    public:
    CongruenceMethod(int a,int c,int m){
        this->a=a,this->c=c,this->m=m;
        r=1;//r0=1とする
        x=(double)r/m;//x0の計算
    }
    //min<x<maxの範囲の疑似乱数xを得る
    double sampling(int min,int max){
        this->r=(a*r+c)%m;
        this->x=(double)r/m;
        return min+x*(max-min);
    }
};

double boxmuller(double mu, double sigma, CongruenceMethod *uniform)
{
    double num1, num2, x1, x2, z1, z2;
    // 0~1の一様乱数生成
    num1 = uniform->sampling(0,1);
    num2 = uniform->sampling(0, 1);
    // ボックスミュラー法
    x1 = sqrt(-2.0 * log(num1)) * cos(2 * pi * num2);
    x2 = sqrt(-2.0 * log(num1)) * sin(2 * pi * num2);
    // 線形変換
    z1 = mu + sigma * x1;
    z2 = mu + sigma * x2;
    //cout << z1 << "," << z2 << endl;
    return z1;
}
double beta(double a,double b,CongruenceMethod *uniform){
    if(a<=0 || b<=0) return -1;
    double alpha=a+b;
    double beta=sqrt((alpha-2)/(2*a*b-alpha));
    if(min(a,b)<=1) beta=max(1/a,1/b);
    double gamma=a+1/beta;
    int loop=1;
    double W;
    while(loop){
        double U1=uniform->sampling(0,1);
        double U2 = uniform->sampling(0, 1);
        double V=beta*log(U1/(1-U1));
        W=a*exp(V);
        loop = (alpha * log(a / (b + W)) + gamma * V - log(4) )< log(U1 * U1 * U2);
    }
    return W/(b+W);
}
int poisson(int lambda,CongruenceMethod *uniform){
    double x=uniform->sampling(0,1);
    double y=exp((double)lambda)*x;
    int n=0;
    while(y>1){
        x=uniform->sampling(0,1);
        y*=x;
        n++;
    }
    return n;
}

int main(){
    //ファイルの作成
    ofstream file1("uniform.csv");
    ofstream file2("norm.csv");
    ofstream file3("beta.csv");
    ofstream file4("poisson.csv");
    int n = 10000;//サンプリング数の指定
    //一様分布の疑似乱数
    int min = 0, max = 10; //一様乱数の範囲
    CongruenceMethod uniform(1229,351570,1664503);
    for(int i=0;i<n;i++){
        file1<<uniform.sampling(min,max)<<"\n";
    }
    //正規分布
    double mu=0,sigma=5;
    for (int i = 0; i < n; i++)
    {
        file2 << boxmuller(mu, sigma, &uniform)<< "\n";
    }
    //ベータ分布
    double a = 2, b = 4;
    for (int i = 0; i < n; i++)
    {
        file3 << beta(a, b, &uniform) << "\n";
    }
    //ポアソン分布
    int lambda=5;
    for (int i = 0; i < n; i++)
    {
        file4 << poisson(lambda, &uniform) << "\n";
    }
    return 0;
}