#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

ZZ log(ZZ a,ZZ g,ZZ p)
{
    ZZ i,j;
    for(i=ZZ(1) ; i<p ; i++)
        if(PowerMod(g,i,p)==a)
            return i;
}

ZZ P_H(ZZ a,ZZ g,ZZ p,Mat<ZZ> &P_1,int l)
{
    ZZ z;
    ZZ Y=a;
    ZZ G=g;
    ZZ O=p-1;
    ZZ P=ZZ(1);
    ZZ X=ZZ(0);
    for(int i=0 ; i<l ; i++)
    {
        for(int j=0 ; j<P_1[i][1] ; j++)
        {
            O/=P_1[i][0];
            z=log(PowerMod(Y,O,p),PowerMod(G,O,p),p);
            Y=(Y*PowerMod(InvMod(G,p),z,p))%p;
            G=PowerMod(G,P_1[i][0],p);
            X=X+P*z;
            P=P*P_1[i][0];
        }
    }
    return X%(p-1);
}

int main()
{
    int l;
    ZZ a,g,p,x;
    cin >> l;
    Mat<ZZ> P;
    P.SetDims(l,2);
    for(int i=0 ; i<l ;i++)
        cin >> P[i][0] >> P[i][1];  ///p=(p1^k1)*...*(pl^kl)
    cin >> a >> g >> p;         ///log_g(x) mod p = a
    x=P_H(a,g,p,P,l);
    cout << x << endl;
    system("PAUSE");
    return 0;
}
