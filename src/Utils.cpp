#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>


long factorial(long n)
{
    if (n==0)
        return 1;
    return n*factorial(n-1);
}


double Cmn(int m, int n)
{
    int k=0, num=1, den=1;
    while(k<m)
    {
        num *= (n-k);
        den *= (k+1);
    }
    return double (num/den);
}

double Pab_k(double csi, int k, int a, int b)
{
    double resp=0.0, a1=0.0, a2=0.0, a3=0.0;

    if (k==0) 
    {
        resp = 1.0;
    }
    else if (k==1)
    {
        resp = 0.5*(a-b+(a+b+2)*csi);
    } 
    else
    {
        a1 = 2*(k+1)*(k+a+b+1)*(2.0*k+a+b);
        a2 = (2*k+a+b+1)*(pow(a,2) - pow(b,2)) + csi*factorial(2*k+a+b+2)/factorial(2*k+a+b-1);
        a3 = 2*(k+a)*(k+b)*(2*k+a+b+2);
        
        resp = (a2*Pab_k(csi,k-1,a,b) - a3*Pab_k(csi,k-2,a,b))/a1;
    }
    return resp;
}

double dPab_k(double csi, int k, int a, int b, int m)
{
    double resp=0.0;
    if (m<=k) 
    {
    resp = (pow(2,-m))*(factorial(k+a+b+m)/factorial(k+a+b))*Pab_k(csi,k-m,a+m,b+m);
    }
    else
    {
        resp = 0.0;
    }
    return resp;
}

double Sfunc(double x2, std::vector<double> resp, int n)
{
    double resp2 = 0.0;
    for (auto i=1; i<=n; i++)
    {
        resp2 = resp2 + pow((x2 - resp[i]),-1);
    }
    return resp2;
}

std::vector<double> newton_raphson(int a, int b, int m)
{
    std::vector<double> resp(m, 0.0);
    double eps = 1E-6, P, dP, S, P_dP;
    double delta = 1.0;
    
    resp[0] =-cos(M_PI*(2*0+1)/2*m);
    for (auto i=0; i<m; i++)
    {
        x2 = -cos(M_PI*(2*i+1)/2*m);
        while (abs(delta) > eps)
        {   
            P  = Pab_k(x2,m,a,b);
            dP = dPab_k(x2,m,a,b,1);
            S  = Sfunc(x2,resp,i);
            P_dP = P/(dP - P*S);

            delta = -P_dP;

            x2 = x2 + delta;
        }
        resp[i] = x2;
    }
    return resp;
}
    

// Jacobi polynomials
double jacobipol(double csi, int order, int a, int b)
{
    std::vector<double > Pvec(order+1, 0.0);
    double a1=0.0,a2=0.0,a3=0.0,a4=0.0;
    
    Pvec[0] = 1.0;
    Pvec[1] = 0.5*(a-b+(a+b+2)*csi);
    for (auto n=1;n<=order-1;n++)
    {
        a1 = 2*(n+1)*(n+a+b+1)*(2*n+a+b);
        a2 = (2*n+a+b+1)*(a^2-b^2);
        a3 = (2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
        a4 = 2*(n+a)*(n+b)*(2*n+a+b+2);
        Pvec[n+1] = ((a2+csi*a3)*Pvec[n]-Pvec[n-1]*a4)/a1;
    }
    return Pvec[order];
}

// Derivative of Jacobi polynomials
double djacobipol(double csi, int order, int a, int b)
{
    if (order == 0)
    {
        return 0.0;
    }
    else
    {
        return 0.5*(a+b+order+1)*jacobipol(csi,order-1,a+1,b+1);
    }
}

// Gauss-Lobatto-Legendre quadrature nodes and weights ---------------------
void NWquadratureGLL(int Qinput, std::vector<double>& nodes, std::vector<double>& weights)
{
    int Q = Qinput-2; // increase for accurate quadrature

    //nodes = zeros(1,Q);
    //std::vector<double> nodes(0.0, Q+2);
    nodes[0] = -1.0;
    nodes[Q+1] = 1.0;


    double tol = 1E-16; // tolerance for root evaluation
    double delta, r, s, jp, djp;

    for (auto k=1; k<=Q; k++)
    {
        delta = tol + 1;
        r = -cos((k-0.5)*M_PI/Q); // initial guess
        if (k > 1)
            r = (r+nodes[k-1])/2; // better initial guess
        
        while (abs(delta) > tol)
        {
            s = 0;
            if (k > 1)
                for (auto i=1; i<=k-1; i++)
                    s = s + 1/(r-nodes[i]);

            jp = jacobipol(r,Q,1,1);
            djp = djacobipol(r,Q,1,1);
            delta = -jp/(djp-s*jp);
            r = r + delta;
        }
        nodes[k] = r; // inner quadrature nodes
    }
    Q = Q + 2;

    double L=0.0;
    for (auto k=0; k<Q; k++)
    {
        L = jacobipol(nodes[k],Q-1,0,0);
        weights[k] = 2/(Q*(Q-1)*pow(L,2));
    }
    // set up output of this routine
    //NW = [nodes, weights];
}

// Lagrange polynomials
// Returns the value of the Lagrange polynomial mode of given INDEX at KSI,
// each polynomial being of order Q, given Q+1 nodes within [-1,1].
// Note that for the nodes and modes, INDEX runs from 0 to Q.
double lagrangepol(double csi, int index, std::vector<double> nodes)
{
    size_t Q = nodes.size() - 1;

    double num=1.0, den=1.0;
    for (auto i=0; i<=Q; i++)
    {
        if (i != index)
        {
            num = num*(csi - nodes[i]);
            den = den*(nodes[index] - nodes[i]);
        }
    }
    return num/den;
}

// Derivative of Lagrange pol.
double dlagrangepol(double csi, int index, std::vector<double> nodes)
{
    size_t Q = nodes.size() - 1;

    // % if (sum(ksi == nodes) == 1) % ksi is one of the nodes
    // %     ksindex = (ksi == nodes)*(0:Q);
    // %     if (ksindex == index)
    // %         num = 0; den = 1;
    // %         for i = 0:Q
    // %             if (i ~= index)
    // %                 num = num + (ksi - nodes(i+1))^-1;
    // %             end
    // %         end
    // %     else
    // %         num = 1; den = 1;
    // %         for i = 0:Q
    // %             if (i ~= index)
    // %                 if (i ~= ksindex)
    // %                     num = num*(ksi - nodes(i+1));
    // %                 end
    // %                 den = den*(nodes(index+1) - nodes(i+1));
    // %             end
    // %         end
    // %     end
    // %else % ksi is not one of the nodes
    double num=0.0, den=1.0, prod=0.0;

    for (auto i=0; i<=Q; i++)
    {
        prod = 0.0;
        if (i != index)
        {
            prod = 1;
            for (auto j=0;j<=Q;j++)
            {
                if (j != index && j != i)
                    prod = prod*(csi - nodes[j]);
            }
            den = den*(nodes[index] - nodes[i]);
        }
        num = num + prod;
    }
    //%end
    return num/den;
}
