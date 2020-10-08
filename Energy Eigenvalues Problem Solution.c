#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int lmin=0;           //Lower boundary of l for energy eigenvalues calculation
int lmax=5;           //Upper boundary of l for energy eigenvalues calculation
double Ei = -30.0;    //Bottom energy to run the experiment
double Ef = -0.00;     //Energy final to end the experiment
int A=244;            //Mass number of the nucleus
int Z=94;             //Atomic number of the nucleus
double h = 0.01;      //Runge kutta step (smaller-->(better accurancy,slower))(default 0.01)
double Estep = 0.001;  //Energy step(smaller->slower, more accurate)

double f(double x, double y, double E, int l, int Vls);
double f1(double x, double y, double yt, double E, int l, int Vls);
double f2(double x, double y, double yt, double E, int l, int Vls);
double f3(double x, double y, double yt, double E, int l, int Vls);
double f4(double x, double y, double yt, double E, int l, int Vls);
double R(int A);

int main()
{
    FILE *fp;

    double E = Ei;      //Eigen value of energy of the nucleus
    int l = lmin;       //Angular momentum of the nucleus
    int Vls;            //Spit-orbit Potential;

    int flag = 0;       //Because a region of energies of the same solution may be found, flag is used to extract the mean value of the solution

    double EF, EI;      //Energy values to calculate mean solution, if continuous solution are found in E+dE
    
    double xe0 ,ye0, yte0; //Initial  Conditions for Runge-Kutta 4
    double xei, yei, ytei; //
    double pya,pyb;        //
    
    fp = fopen("RK4.dat", "w");
        
    while(l<=lmax){

	printf("Energy eigenvalues for element with A=%d and Z=%d\n", A, Z);
	printf("l = %d:\n",l);

	for(Vls=0; Vls<2; Vls++){
	    
	    while(E<Ef){

		xe0 = 0.00000000000000001;        //initial value of x0(=0.00000000001 for nan error, computational reasons)
		ye0 = 0.0;                  //initial value of y(x0)(=0.0)
		yte0 = 0.1;                   //initial value of y'(x0)(can be whaterver as it only changes the amplitude of the y(x). Small for small amplitude)

		while((xe0)<R(A))           //Runge kutta until R(A)
		    {
			xei = xe0 + h;
	    
			ytei = yte0 + 0.166666666*h*(f1(xe0, ye0, yte0, E, l, Vls)+2.0*f2(xe0, ye0, yte0, E, l, Vls)+2.0*f3(xe0, ye0, yte0, E, l, Vls)+f4(xe0, ye0, yte0, E, l, Vls));

			yte0 = ytei;

			yei = ye0 + yte0*h + 0.166666666*h*h*(f1(xe0, ye0, yte0, E, l, Vls)+f2(xe0, ye0, yte0, E, l, Vls)+f3(xe0, ye0, yte0, E, l, Vls));
	    
			xe0 = xei;
			pya = ye0 = yei;
		    }

		{                             //Runge Kutta step to R(A)+h
		    xei = xe0 + h;
	
		    ytei = yte0 + 0.166666666*h*(f1(xe0, ye0, yte0, E, l, Vls)+2.0*f2(xe0, ye0, yte0, E, l, Vls)+2.0*f3(xe0, ye0, yte0, E, l, Vls)+f4(xe0, ye0, yte0, E, l, Vls));
	
		    yei = ye0 + yte0*h + 0.166666666*h*h*(f1(xe0, ye0, yte0, E, l, Vls)+f2(xe0, ye0, yte0, E, l, Vls)+f3(xe0, ye0, yte0, E, l, Vls));
	
		    xe0 = xei;
		    pyb = ye0 = yei;
		    yte0 = ytei;

		}
	
		if(pya*pyb < 0.0)              //Is y( R(A) ) a root (Newton's method (Also known to us, as bolzano theorem))
		    {
			if(flag == 0) //Are we still in E+dE, or a new solution?
			    {
				flag = 1;
				EI = E;
			    }
			
		    }
		else
		    {
			if(flag == 1) //Is this the end of E+dE?
			    {
				flag = 0;
				EF = E - Estep;
				printf("%lf\n", (EF+EI)/2.0);
				fprintf(fp,"%lf\n", (EF+EI)/2.0);
			    }
		    }
		    
		E = E + Estep;                  //Increase energy and try again(ENERGY STEP = 0.01(default))
	    }
	    flag = 0;
	    E = Ei;
	}
	fprintf(fp, "\n");
	l = l + 1;                         //Increase angular momentum
    }

    fclose(fp);
    
    return 0;
}


double f(double x, double y, double E, int l, int Vls)                      // du'/dx=f(x,y), u'=du/dx (The nuclear potential of the nucleus)
{
    if(Vls == 0 && l != 0)
	{
	    return ((-0.2*((double)l)/( exp((x-R(A))/0.6) + 1.0 ) - 43.0/( exp((x-R(A))/0.6) + 1.0 ) + 20.0*((double)l*((double)l+1.0))/(x*x) + 1.44*((double)Z-1.0)/x)-E)*y;
	}
    if(Vls == 1 && l != 0)
	{
	    return ((0.2*((double)l+1.0)/( exp((x-R(A))/0.6) + 1.0 ) - 43.0/( exp((x-R(A))/0.6) + 1.0 ) + 20.0*((double)l*((double)l+1.0))/(x*x) + 1.44*((double)Z-1.0)/x)-E)*y;
	}
    if(Vls == 0 && l == 0)
	{
	    return ((-43.0/( exp((x-R(A))/0.6) + 1.0 ) + 20.0*((double)l*((double)l+1.0))/(x*x) + 1.44*((double)Z-1.0)/x)-E)*y;
	}
    if(Vls == 1 && l == 0)
	{
	    return 0;
	}
    exit(-1);
}
double f1(double x, double y, double yt, double E, int l, int Vls)          //Runge kutta f1
{
    return f(x, y, E, l, Vls);
}
double f2(double x, double y, double yt, double E, int l, int Vls)          //Runge kutta f2
{
    return f(x+0.5*h, y+0.5*h*yt, E, l, Vls);
}
double f3(double x, double y, double yt, double E, int l, int Vls)          //Runge kutta f3
{
    return f(x+0.5*h, y+yt*h*0.5+f1(x, y, yt, E, l, Vls)*0.25*h*h, E, l, Vls);
}
double f4(double x, double y, double yt, double E, int l, int Vls)          //Runge kutta f4
{
    return f(x+h, y+yt*h+f2(x, y, yt, E, l, Vls)*0.5*h, E, l, Vls);
}
double R(int A)                                             //Nuclear radius function
{
    return 1.2*pow(A, 0.333);
}
