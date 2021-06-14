#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define Maxpnt  5000 /* the maximum number of points. this can be changed for a larger system */
#define dim 3
#define G 6.67408e-11 /* units of m3 /kg /s2 */
#define kpc 3.086e19 /*conversion of meters to kiloparsecs */
#define Msolar 1.98855e30 /* mass of the sun in kg */
#define Pi 3.14159265358979323846
#define bmass 5.683e26 /* the mass of the central  body */
#define lmass 1e24 /* the mass of the orbiting bodies*/
/*#define baxis sqrt(1-(ecc*ecc))*axis */
#define ecc sqrt(1-(baxis*baxis)/(axis*axis))

void main(argc, argv)
int argc;
char *argv[];
{
    int i,sweep,burn,meas,ind,u=0,j,N;
    double pos[Maxpnt],x[Maxpnt],y[Maxpnt],z[Maxpnt],dr, dx,dy,dz, pos1,x1,y1,z1,r,delta,ds,l=0,acc=0,NN,x2,z2,y2,pos2,e=0,e2=0;
    srand((unsigned)time(NULL));

    double h=1.0;
    double m=1.0;
    double w=1.0;
    double tau=1000;
    N = 2000;
    NN= 2000.0;
    double a=tau/NN;
    delta=.6;
    double c=pow(a*w,2.0);
    double w2=pow(w,2.0);
    double ee=0.0,ee2=0.0;

    burn=100000;
    meas=1000000;
    ind=1000;

    FILE *poshist;
    FILE *pos2hist;
    poshist = fopen("./poshist.c","w");
    pos2hist = fopen("./pos2hist.c","w");


    for(i=0;i<N;i++)
    {
        x[i]=0.0;
        y[i]=0.0;
        z[i]=0.0;
    }
    x[N]=x[0];
    y[N]=y[0];
    z[N]=z[0];

    x[-1]=x[N-1];
    y[-1]=y[N-1];
    z[-1]=z[N-1];


    for(sweep=0;sweep<meas;sweep++)
    {
        for(i=0;i<N;i++)
        {
            dx=delta*2*(((double)rand()/(double)RAND_MAX)-.5);
            dy=delta*2*(((double)rand()/(double)RAND_MAX)-.5);
            dz=delta*2*(((double)rand()/(double)RAND_MAX)-.5);

            x1=x[i]+dx;
            y1=y[i]+dy;
            z1=z[i]+dz;

            //ds= dx *(2*u[i] + dx + g*(u[i-1] + u[i+1]));


            ds= (dx*(dx*(2*m + c) + 2*c*x[i] - 2*m*(x[i-1] - 2*x[i] + x[i+1])))/(2*a);
            ds+= (dy*(dy*(2*m + c) + 2*c*y[i] - 2*m*(y[i-1] - 2*y[i] + y[i+1])))/(2*a);
            ds+=(dz*(dz*(2*m + c) + 2*c*z[i] - 2*m*(z[i-1] - 2*z[i] + z[i+1])))/(2*a);


            //ds= (m/(2*a*a))* (2*dx*dx + 2*dy*dy + 2*dz*dz - 2*dx*x[i-1] + 4*dx*x[i] - 2*dx*x[i+1] - 2 dy y2 + 4 dy y3 - 2 dy y4 - 2 dz z2 + 4 dz z3 + a^2 w^2 Sqrt[x3^2 + y3^2 + z3^2] - a^2 w^2 Sqrt[(dx + x3)^2 + (dy + y3)^2 + (dz + z3)^2] - 2 dz z4)

            //printf("ds= %e \n",exp(-ds));

            if(ds<0){
                x[i]=x1;
                y[i]=y1;
                z[i]=z1;
                acc+=1;
            }
            else{
                if( ((double)rand()/(double)RAND_MAX) <= exp(-ds) )
                {
                    x[i]=x1;
                    y[i]=y1;
                    z[i]=z1;
                    acc+=1;
                }
            }

            l+=1;
            r=acc/l;
            x[N]=x[0];
            y[N]=y[0];
            z[N]=z[0];

            x[-1]=x[N-1];
            y[-1]=y[N-1];
            z[-1]=z[N-1];
        }

        if(sweep == burn + u*ind){
            x2=0.0;
            pos1=0.0;
            pos2=0.0;
            e=0;e2=0;
            for(j=0;j<N;j++){
                x2=x[j]*x[j];
                pos2+=pow(x[j],2.0);
                pos2+=pow(y[j],2.0);
                pos2+=pow(z[j],2.0);
                pos1+=x[j];
                pos1+=y[j];
                pos1+=z[j];
                fprintf(poshist,"%f  %f  %f \n",x[j], y[j], z[j]);
                //fprintf(pos2hist,"%f \n",u2);
                e+=m*w2*x2;
                e2+=pow(m*w2*x2,2.0);
            }
            pos1=pos1/NN;
            pos2=pos2/NN;
            ee+=m*pow(w,2.0)*pos2;
            ee2+=pow(m*pow(w,2.0)*pos2,2.0);
            //fprintf(poshist,"%f \n",pos1);
            fprintf(pos2hist,"%f \n",pos2);
            u+=1;
            printf("%i   acceptance ratio = %f\n",u,r );
        }

    }
    ee=ee/u;
    ee2=ee2/u;
    double er=sqrt((ee2-pow(ee,2.0)))/sqrt(sqrt(u));
    printf("the ground state energy is %f \nthe error in the energy is %f ",ee,er);

}
