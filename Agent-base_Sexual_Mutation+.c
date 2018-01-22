/*w
 Author Shota Shibasaki
 Last update 15/01/2018
 Agent-based Sexual model
 without mutation
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include"MT.h"
/*
 necessary for Mersenne Twister. Copy rght (C) 1997-2002 Makoto Matsumoto and Tukuji NIshimura.
 http://www.sat.t.u-tokyo.ac.jp/~omi/code/MT.h
 */

double *malloc_vector(int n);
int *malloc_vectori(int n);
void free_vector(double *a);
void free_vectori(int *a);
double **malloc_matrix(int nrow, int ncol);
int **malloc_matrixi(int nrow, int ncol);
void free_matrix(double **a, int nrow);
void free_matrixi(int **a, int nrow);
char fname[256];
char sname1[256];
char sname2[256];

/*main function*/
int main(void)
{
    clock_t time_st, time_end;
    int i; int j; int k; int l; int I; int J;
    int strain;
    strain=2;//Two Strain
    int totstrain=6;//total number of strain
    double rand;
    int the;
    int rep;
    int c; int r;
    int t=0; int T=100000;//simulation from 0 to T.
    int Nt=300;
    int ntot=Nt*Nt;
    int th=1;
    double theta;//probability of macrocyst formation between the mating types
    double prob;
    double max;
    int maxcount;
    int *maxi;
    maxi=malloc_vectori(9);
    int *maxj;
    maxj=malloc_vectori(9);
    int *count;//number of colonies of strategy i
    count=malloc_vectori(totstrain);
    for (i=1;i<=totstrain;i++)
    {
        count[i]=0;
    }
    int tr; int Tr=10;//number of trial
    int **n;// strain type
    n = malloc_matrixi(Nt,Nt);
    int **m;// strain type used in update
    m=malloc_matrixi(Nt,Nt);
    int **x;//strategy
    x=malloc_matrixi(Nt,Nt);
    int **s;//mating type
    s=malloc_matrixi(Nt,Nt);
    double **fit;//fitness
    fit=malloc_matrix(Nt,Nt);
    double **phi;//fitness used in update
    phi=malloc_matrix(Nt,Nt);
    int **idlabel;
    idlabel=malloc_matrixi(Nt,Nt);
    int *idc;//individual collum
    idc=malloc_vectori(ntot);
    int *idr;//individual row
    idr=malloc_vectori(ntot);
    int idcount;
    int opplabel;
    double **benefit;
    benefit=malloc_matrix(ntot,8);//benefit of each individuals
    int opponent;//number of opponent
    // set the payoff matrices
    double e0=0.000001;
    double **A; A=malloc_matrix(3,3);//matrix A, fruiting body foramtoin
    double hat=0.5;
    double a1=1.0; double b1=2.7; double e1=0.00001;
    int b1_count;
    A[1][1]=a1;
    A[1][2]=hat*a1;
    A[1][3]=0.0;
    A[2][1]=e0;
    A[2][2]=e0;
    A[2][3]=e0;
    A[3][1]=b1;
    A[3][2]=hat*e1;
    A[3][3]=e1;
    double **B; B=malloc_matrix(3,3);//matrix B, macrocyst formation
    double a2=0.6; double b2=a2+0.6; double e2=0.00001;
    int a2_counter;
    B[1][1]=a2;
    B[1][2]=0.0;
    B[1][3]=a2;
    B[2][1]=b2;
    B[2][2]=e2;
    B[2][3]=b2;
    B[3][1]=a2;
    B[3][2]=0.0;
    B[3][3]=a2;
    double **P;//probability of update
    P=malloc_matrix(ntot,9);
    double K=0.1;//uncertainty of the adoptation
    //K=0.0 means  update with the highest strategy but Nan can be produced in C
    double NP;//normirizing P
    double L;
    //using update process
    int *oppi; int *oppj;
    oppi=malloc_vectori(9);
    oppj=malloc_vectori(9);
    // end the matrices setting
    
   // printf("Debug 0\n");
    time_st=clock();
    theta=0.8;
    printf("  theta = %.1f\n",theta);
    for (tr=1;tr<=Tr;tr++)//trial loop
    {
        t=0;
        //random seed
        init_genrand((unsigned)time(NULL)+tr);
        /*set initial condition*/
        for (i=1;i<=totstrain;i++)
        {
            count[i]=0;
        }
        idcount=1;
        for (r=1;r<=Nt;r++)
        {
            for (c=1;c<=Nt;c++)
            {
                //set the strategy
                n[c][r]=genrand_int32()%(2-1+1)+1;//random number from 1 to strain
                //only C1 and C2 at first
                    //random int from A to B => genrand_int32()%(B-A+1)+A
                if (n[c][r]==2)
                {
                    n[c][r]=4;//C2
                }
                //strategy
                if (n[c][r]<4)
                {
                    s[c][r]=1;//mating type1
                    x[c][r]=n[c][r];//strategy
                }
                else
                {
                    s[c][r]=2;//mating type 2
                    x[c][r]=n[c][r]-3;//strategy
                }
                m[c][r]=n[c][r];
                count[n[c][r]]+=1;
                fit[c][r]=0;
                idc[idcount]=c;
                idr[idcount]=r;
                idlabel[c][r]=idcount;
                idcount=idcount+1;
            }
        }
        // printf("Debug 1\n");
        //save initial condition
        FILE *fp;
        sprintf(fname,"mix_+mutation_theta-%.1f_K%.1f_trial%d.csv",theta,K,tr);
        fp=fopen(fname,"w");
        fprintf(fp,"%d,",t);
        for (i=1;i<=totstrain;i++)
        {
            if(i==totstrain)
            {
                fprintf(fp,"%d\n",count[i]);
            }
            else
            {
                fprintf(fp,"%d,",count[i]);
            }
        }
        fclose(fp);
        //printf("debug2\n");
        
        //start dynamics
        for (t=1;t<=T;t++)//time loop
        {
            //reset the benefit
            for (i=1;i<=ntot;i++)
            {
                for (j=1;j<=8;j++)
                {
                    benefit[i][j]=-1;
                }
            }
            //interacting with neighbors
            for (idcount=1;idcount<=ntot;idcount++)
            {
                r=idr[idcount];
                c=idc[idcount];
                //printf("player %d, position (%d, %d)\n",idcount,c,r);
                opponent=1;
                for (i=c-1;i<=c+1;i++)
                {
                    if(i<1)
                    {
                        I=Nt;
                    }
                    else if (i>Nt)
                    {
                        I=1;
                    }
                    else
                    {
                        I=i;
                    }
                    for(j=r-1;j<=r+1;j++)
                    {
                        if(j<1)
                        {
                            J=Nt;
                        }
                        else if (j>Nt)
                        {
                            J=1;
                        }
                        else
                        {
                            J=j;
                        }
                        //  printf("vs opponent %d position (%d, %d)\n",opponent, I,J);
                        opplabel=idlabel[I][J];
                        if (idcount!=opplabel)//interaction does not occuer with inself
                        {
                            if (benefit[idcount][opponent]<0)//yet interacted
                            {
                                if(s[c][r]==s[I][J])//within mating types game
                                {
                                    benefit[idcount][opponent]=A[x[c][r]][x[I][J]];
                                    benefit[opplabel][9-opponent]=A[x[I][J]][x[c][r]];
                                }
                                else//between mating types game
                                {
                                    prob=genrand_real1();
                                    if (prob<=theta)
                                    {
                                        benefit[idcount][opponent]=B[x[c][r]][x[I][J]];
                                        benefit[opplabel][9-opponent]=B[x[I][J]][x[c][r]];
                                    }
                                    else
                                    {
                                        benefit[idcount][opponent]=A[x[c][r]][x[I][J]];
                                        benefit[opplabel][9-opponent]=A[x[I][J]][x[c][r]];
                                    }
                                }
                                // printf("Payoff %.1f\n",benefit[idcount][opponent]);
                                fit[c][r]+=benefit[idcount][opponent];
                                fit[I][J]+=benefit[opplabel][9-opponent];
                                opponent=opponent+1;
                            }
                            else
                            {
                                // printf("ALREADY INTERACTED\n");
                                //printf("Payoff %.1f\n",benefit[idcount][opponent]);
                                opponent=opponent+1;
                            }
                        }
                        else
                        {
                            //printf("itself PASS ME!\n");
                        }
                    }
                }
                //printf("agent %d fitness %.1f\n\n",idcount,fit[c][r]);
                phi[c][r]=fit[c][r];
            }
            // printf("Debug2.5\n");
                
            //update
            for(idcount=1;idcount<=ntot;idcount++)
            {
                //printf("player %d",idcount);
                r=idr[idcount];
                c=idc[idcount];
                //printf("player (%d, %d)\n",c,r);
                opponent=1;
                NP=0;
                for(i=c-1;i<=c+1;i++)
                {
                    if(i<1)
                    {
                        I=Nt;
                    }
                    else if (i>Nt)
                    {
                        I=1;
                    }
                    else
                    {
                        I=i;
                    }
                    for (j=r-1;j<=r+1;j++)
                    {
                        if(j<1)
                        {
                            J=Nt;
                        }
                        else if (j>Nt)
                        {
                            J=1;
                        }
                        else
                        {
                            J=j;
                        }
                        oppi[opponent]=I;
                        oppj[opponent]=J;
                        //printf("vs (%d, %d)\n",oppi[opponent],oppj[opponent]);
                        P[idcount][opponent]=exp((fit[oppi[opponent]][oppj[opponent]]-fit[c][r])/K);
                        NP=NP+P[idcount][opponent];
                        opponent=opponent+1;
                    }
                }
                if (isnan(NP))
                {
                    printf("Nan ! Error in sum of exp \n");
                }
                //printf("end cal prob\n");
                for (k=1;k<=9;k++)
                {
                    P[idcount][k]=P[idcount][k]/NP;
                    if (isnan(P[idcount][k]))
                    {
                        printf("Nan ! Error in softmax function\n");
                        return 1;
                    }
                    else
                    {
                        //printf("prob %.3f\n",P[idcount][k]);
                    }
                }
                //update to the strategy of one neighbor
                prob=genrand_real1();
                //printf("probability %.6f\n",prob);
                L=0.0;
                for (l=1;l<=9;l++)
                {
                    //printf("fitness %.6f \n",P[idcount][l]);
                    if (prob<=L+P[idcount][l])
                    {
                        //update with label i
                        I=oppi[l];
                        J=oppj[l];
                        //printf("Update from %d to (%d, %d)\n",idcount,I,J);
                        m[c][r]=n[I][J];
                        //phi[c][r]=fit[oppi[i]][oppj[i]];
                        L=0.0;
                        break;
                    }
                    else
                    {
                        L=L+P[idcount][l];
                    }
                }
            }// end update
                
            //printf("Debug 3\n");
            // strategy update
            for (i=1;i<=totstrain;i++)
            {
                count[i]=0;
            }
            for (r=1;r<=Nt;r++)
            {
                for(c=1;c<=Nt;c++)
                {
                    n[c][r]=m[c][r];
                    if(n[c][r]<4)
                    {
                        x[c][r]=n[c][r];
                        s[c][r]=1;
                    }
                    else
                    {
                        x[c][r]=n[c][r]-3;
                        s[c][r]=2;
                    }
                    fit[c][r]=0;
                    count[n[c][r]]+=1;
                }
            }//end updating strategy and count
            //printf("Debug 4\n");
            if (t%500==0)//mutation occurs in each 500 time steps
            {
                //choose an agent randomly for mutation
                i=genrand_int32()%(Nt-1+1)+1;
                j=genrand_int32()%(Nt-1+1)+1;
                l=n[i][j];
                while (n[i][j]==l)
                {
                    l=genrand_int32()%(totstrain-1+1)+1;//randomly updated
                }
                //printf("mutation at (%d, %d). from %d to %d\n",i,j,n[i][j],l);
                count[n[i][j]]-=1;
                count[l]+=1;
                n[i][j]=l;
            }
            if(t%100==0)//save loop
            {
                //save the strategies
                FILE *fp;
                sprintf(fname,"mix_+mutation_theta-%.1f_K%.1f_trial%d.csv",theta,K,tr);
                fp=fopen(fname,"a");
                fprintf(fp,"%d,",t);
                for (i=1;i<=totstrain;i++)
                {
                    if(i==totstrain)
                    {
                        fprintf(fp,"%d\n",count[i]);
                    }
                    else
                    {
                        fprintf(fp,"%d,",count[i]);
                    }
                }
                fclose(fp);
            }//save loop
            //mutaion
        }//time loop
    }//trial loop
    printf("end the simulation\n");
    return 0;
}
/*Dynamic storage alloaction*/
/*vecotr ver*/
double *malloc_vector(int n)
{
    double *a;
    if ((a=malloc(sizeof(double)*n))==NULL)
    {
        printf("unable to keep memory for vector, sorry\n");
        exit(1);
    }
    return (a-1);/*row number stratswith number 1*/
}
void free_vector(double *a)
/*free the memory*/
{
    free(a+1);
}
int *malloc_vectori(int n)
{
    int *a;
    if ((a=malloc(sizeof(double)*n))==NULL)
    {
        printf("unable to keep memory for vector, sorry\n");
        exit(1);
    }
    return (a-1);/*row number stratswith number 1*/
}
void free_vectori(int *a)
/*free the memory*/
{
    free(a+1);
}
/*matrix ver*/
double **malloc_matrix(int nrow, int ncol)
{
    double **a;
    int i;
    if ((a=malloc(nrow*sizeof(double *)))==NULL)
    {
        printf("unable to keep memory for matrix sorry");
        exit(1);
    }
    a = a-1;
    for (i=1;i<=nrow;i++) a[i]=malloc(ncol*sizeof(double));
    for (i=1;i<=nrow;i++) a[i]=a[i]-1;/* move row direction*/
        
        return a;
}
void free_matrix(double **a, int nrow)
{
    int i;
    for ( i=1; i<=nrow; i++) free((void *)(a[i]+1));
    free((void *)(a+1));
}
/*matrix ver*/
int **malloc_matrixi(int nrow, int ncol)
{
    int **a;
    int i;
    if ((a=malloc(nrow*sizeof(double *)))==NULL)
    {
        printf("unable to keep memory for matrix sorry");
        exit(1);
    }
    a = a-1;
    for (i=1;i<=nrow;i++) a[i]=malloc(ncol*sizeof(double));
    for (i=1;i<=nrow;i++) a[i]=a[i]-1;/* move row direction*/
    
    return a;
}
void free_matrixi(int **a, int nrow)
{
    int i;
    for ( i=1; i<=nrow; i++) free((void *)(a[i]+1));
    free((void *)(a+1));
}