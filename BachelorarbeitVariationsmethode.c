#include <stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

int LR(int n, double** A);
void ausgabe_Vektor(double* b, int n);
void ausgabe_Matrix(double** matrix, int n, int m);
int vwsubs(int n, double** B, double* b, double* d);
int rwsubs(int n, double** B, double* b);
double pot(double d, int n);

int main(void)
{

    double** matrixDH;
	double* vektornegH;
	double* loesungsvektor;
	double* kopieloesungsvektor;
	double matrixA_diagonale;
    double matrixA_nichtdiagonale;
    int wiederholungen;//Anzahl der Widerholungen beim Newton-Verfahren
    int k;//Zaehlvariable
	int n;// Anzahl der Zeitintervalle - 1
	double R;//Grenzen des Definitionsbereichs
	double h;//Schrittweite
	double tau;//Schrittweite der Diskretisierung des partiellen Diffgl

	//Einstellung der Parameter
	wiederholungen = 10;
	n = 100;
	R = 5;
	h = (R+R)/(n+1);
	tau = 1.5;
    k = 0;

    //Speicherplatz für Matrizen und Vektoren wird festgelegt
    loesungsvektor = (double*)calloc(n, sizeof(double));
	ausgabe_Vektor(loesungsvektor, n);
    matrixDH = (double**)calloc(n, sizeof(double*));
    for(int i=0;i<n;i++)
    {
        matrixDH[i] = (double*)calloc(3, sizeof(double));
    }
    kopieloesungsvektor = (double*)calloc(n, sizeof(double));
    vektornegH = (double*)calloc(n, sizeof(double));

    //Anfangswert fuer das Newton_Verfahren wird festgelegt
    for(int i=0;i<n;i++)
	{
		loesungsvektor[i] = sin(M_PI/(2*R)*(-R+(i+1)*h))- (-R+(i+1)*h)/R;
	}

    //Berechnung der Matrix A
	matrixA_diagonale = 2/h+(1/tau-(double) 1/2)*2*h/3;
	matrixA_nichtdiagonale = -1/h+(1/tau-(double) 1/2)*1/6*h;

    //Beginn des Newton-Verfahrens
	while(k<wiederholungen)
    {
        //Kopie des Loesungsvektors wird erstellt
        for(int i=0;i<n;i++)
            {
            kopieloesungsvektor[i] = loesungsvektor[i];
        }

        //DH wird erstellt
        for(int i=0;i<n;i++)
        {
            if(i>0 && i<(n-1)) matrixDH[i][1] = matrixA_diagonale+3/2*(1/(3*(-loesungsvektor[i-1]/h+loesungsvektor[i]/h+1/R))*pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 3)-2/(3*h*(-loesungsvektor[i-1]/h+loesungsvektor[i]/h+1/R))*(1/(4*(-loesungsvektor[i-1]/h+loesungsvektor[i]/h+1/R))*pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 4)-1/(20*h*pot(-loesungsvektor[i-1]/h+loesungsvektor[i]/h+1/R, 2))*(pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 5)-pot(loesungsvektor[i-1]+(-R+i*h)/R, 5)))-1/(3*(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R))*pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 3)+2/(3*h*(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R))*(-1/(4*(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R))*pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 4)+1/(20*h*pot(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R, 2))*(pot(loesungsvektor[i+1]+(-R+(i+2)*h)/R, 5)-pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 5))));
            if(i<(n-1)) matrixDH[i][2] = matrixA_nichtdiagonale-(1/(2*(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R))*(-1/(h*4*(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R))*(pot(loesungsvektor[i+1]+(-R+(i+2)*h)/R, 4)+pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 4))+1/(10*h*h*pot(-loesungsvektor[i]/h+loesungsvektor[i+1]/h+1/R, 2))*(pot(loesungsvektor[i+1]+(-R+(i+2)*h)/R, 5)-pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 5))));
            if(i>0) matrixDH[i][0] = matrixDH[i-1][2];
        }
        matrixDH[0][1] = matrixA_diagonale+3/2*(1/(3*(loesungsvektor[0]/h+1/R))*pot(loesungsvektor[0]+(-R+h)/R, 3)-2/(3*h*(loesungsvektor[0]/h+1/R))*(1/(4*(loesungsvektor[0]/h+1/R))*pot(loesungsvektor[0]+(-R+h)/R, 4)-1/(20*h*pot(loesungsvektor[0]/h+1/R, 2))*(pot(loesungsvektor[0]+(-R+h)/R, 5)+1))-1/(3*(-loesungsvektor[0]/h+loesungsvektor[1]/h+1/R))*pot(loesungsvektor[0]+(-R+h)/R, 3)+2/(3*h*(-loesungsvektor[0]/h+loesungsvektor[1]/h+1/R))*(-1/(4*(-loesungsvektor[0]/h+loesungsvektor[1]/h+1/R))*pot(loesungsvektor[0]+(-R+h)/R, 4)+1/(20*h*pot(-loesungsvektor[0]/h+loesungsvektor[1]/h+1/R, 2))*(pot(loesungsvektor[1]+(-R+2*h)/R, 5)-pot(loesungsvektor[0]+(-R+h)/R, 5))));
        matrixDH[n-1][1] = matrixA_diagonale+3/2*(1/(3*(-loesungsvektor[n-2]/h+loesungsvektor[n-1]/h+1/R))*pot(loesungsvektor[n-1]+(-R+n*h)/R, 3)-2/(3*h*(-loesungsvektor[n-2]/h+loesungsvektor[n-1]/h+1/R))*(1/(4*(-loesungsvektor[n-2]/h+loesungsvektor[n-1]/h+1/R))*pot(loesungsvektor[n-1]+(-R+n*h)/R, 4)-1/(20*h*pot(-loesungsvektor[n-2]/h+loesungsvektor[n-1]/h+1/R, 2))*(pot(loesungsvektor[n-1]+(-R+n*h)/R, 5)-pot(loesungsvektor[n-2]+(-R+(n-1)*h)/R, 5)))-1/(3*(-loesungsvektor[n-1]/h+1/R))*pot(loesungsvektor[n-1]+(-R+n*h)/R, 3)+2/(3*h*(-loesungsvektor[n-1]/h+1/R))*(-1/(4*(-loesungsvektor[n-1]/h+1/R))*pot(loesungsvektor[n-1]+(-R+n*h)/R, 4)+1/(20*h*pot(-loesungsvektor[n-1]/h+1/R, 2))*(1-pot(loesungsvektor[n-1]+(-R+n*h)/R, 5))));
        matrixDH[0][0] = 0;
        matrixDH[n-1][2] = 0;

        //-H wird erstellt
        vektornegH[0] = -(loesungsvektor[0]*matrixA_diagonale+matrixA_nichtdiagonale*loesungsvektor[1]-((2*R)/(tau*h*M_PI)*(-h*cos(M_PI/(2*R)*(-R+h))+(2*R)/M_PI*sin(M_PI/(2*R)*(-R+h))-(2*R)/M_PI*sin(M_PI/(2*R)*(-R)))-1/(2*h*4*(loesungsvektor[0]/h+1/R))*(h*pot(loesungsvektor[0]+(-R+h)/R, 4)-1/(5*(loesungsvektor[0]/h+1/R))*(pot(loesungsvektor[0]+(-R+h)/R, 5)-pot(-R/R, 5)))-1/h*(1/tau-(double) 1/2)*(pot(-R+h, 3)/(3*R)-((-R)*pot(-R+h, 2))/(2*R)+pot(-R, 3)/(6*R)))-(-1*((2*R)/(tau*h*M_PI)*(-h*cos(M_PI/(2*R)*(-R+2*h))+(2*R)/M_PI*sin(M_PI/(2*R)*(-R+2*h))-(2*R)/M_PI*sin(M_PI/(2*R)*(-R+h)))-1/(2*h*4*(loesungsvektor[1]/h-loesungsvektor[0]/h+1/R))*(h*pot(loesungsvektor[1]+(-R+2*h)/R, 4)-1/(5*(loesungsvektor[1]/h-loesungsvektor[0]/h+1/R))*(pot(loesungsvektor[1]+(-R+2*h)/R, 5)-pot(loesungsvektor[0]+(-R+h)/R, 5)))-1/h*(1/tau-(double) 1/2)*(pot(-R+2*h, 3)/(3*R)-((-R+h)*pot(-R+2*h, 2))/(2*R)+pot(-R+h, 3)/(6*R)))+1/tau*(-cos(M_PI/(2*R)*(-R+2*h))*2*R/M_PI+cos(M_PI/(2*R)*(-R+h))*2*R/M_PI)-1/(8*(loesungsvektor[1]/h-loesungsvektor[0]/h+1/R))*(pot(loesungsvektor[1]+(-R+2*h)/R, 4)-pot(loesungsvektor[0]+(-R+h)/R, 4))-(1/tau-(double) 1/2)*(pot(-R+2*h, 2)/(2*R)-pot(-R+h, 2)/(2*R))));
        for(int i=1;i<n-1;i++)
        {
            vektornegH[i] = -(loesungsvektor[i]*matrixA_diagonale+matrixA_nichtdiagonale*(loesungsvektor[i-1]+loesungsvektor[i+1])-((2*R)/(tau*h*M_PI)*(-h*cos(M_PI/(2*R)*(-R+(i+1)*h))+(2*R)/M_PI*sin(M_PI/(2*R)*(-R+(i+1)*h))-(2*R)/M_PI*sin(M_PI/(2*R)*(-R+i*h)))-1/(2*h*4*(loesungsvektor[i]/h-loesungsvektor[i-1]/h+1/R))*(h*pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 4)-1/(5*(loesungsvektor[i]/h-loesungsvektor[i-1]/h+1/R))*(pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 5)-pot(loesungsvektor[i-1]+(-R+i*h)/R, 5)))-1/h*(1/tau-(double) 1/2)*(pot(-R+(i+1)*h, 3)/(3*R)-((-R+i*h)*pot(-R+(i+1)*h, 2))/(2*R)+pot(-R+i*h, 3)/(6*R)))-(-1*((2*R)/(tau*h*M_PI)*(-h*cos(M_PI/(2*R)*(-R+(i+2)*h))+(2*R)/M_PI*sin(M_PI/(2*R)*(-R+(i+2)*h))-(2*R)/M_PI*sin(M_PI/(2*R)*(-R+(i+1)*h)))-1/(2*h*4*(loesungsvektor[i+1]/h-loesungsvektor[i]/h+1/R))*(h*pot(loesungsvektor[i+1]+(-R+(i+2)*h)/R, 4)-1/(5*(loesungsvektor[i+1]/h-loesungsvektor[i]/h+1/R))*(pot(loesungsvektor[i+1]+(-R+(i+2)*h)/R, 5)-pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 5)))-1/h*(1/tau-(double) 1/2)*(pot(-R+(i+2)*h, 3)/(3*R)-((-R+(i+1)*h)*pot(-R+(i+2)*h, 2))/(2*R)+pot(-R+(i+1)*h, 3)/(6*R)))+1/tau*(-cos(M_PI/(2*R)*(-R+(i+2)*h))*2*R/M_PI+cos(M_PI/(2*R)*(-R+(i+1)*h))*2*R/M_PI)-1/(8*(loesungsvektor[i+1]/h-loesungsvektor[i]/h+1/R))*(pot(loesungsvektor[i+1]+(-R+(i+2)*h)/R, 4)-pot(loesungsvektor[i]+(-R+(i+1)*h)/R, 4))-(1/tau-(double) 1/2)*(pot(-R+(i+2)*h, 2)/(2*R)-pot(-R+(i+1)*h, 2)/(2*R))));
        }
        vektornegH[n-1] = -(loesungsvektor[n-1]*matrixA_diagonale+matrixA_nichtdiagonale*loesungsvektor[n-2]-((2*R)/(tau*h*M_PI)*(-h*cos(M_PI/(2*R)*(-R+n*h))+(2*R)/M_PI*sin(M_PI/(2*R)*(-R+n*h))-(2*R)/M_PI*sin(M_PI/(2*R)*(-R+(n-1)*h)))-1/(2*h*4*(loesungsvektor[n-1]/h-loesungsvektor[n-2]/h+1/R))*(h*pot(loesungsvektor[n-1]+(-R+n*h)/R, 4)-1/(5*(loesungsvektor[n-1]/h-loesungsvektor[n-2]/h+1/R))*(pot(loesungsvektor[n-1]+(-R+n*h)/R, 5)-pot(loesungsvektor[n-2]+(-R+(n-1)*h)/R, 5)))-1/h*(1/tau-(double) 1/2)*(pot(-R+n*h, 3)/(3*R)-((-R+(n-1)*h)*pot(-R+n*h, 2))/(2*R)+pot(-R+(n-1)*h, 3)/(6*R)))-(-1*((2*R)/(tau*h*M_PI)*(-h*cos(M_PI/(2*R)*(-R+(n+1)*h))+(2*R)/M_PI*sin(M_PI/(2*R)*(-R+(n+1)*h))-(2*R)/M_PI*sin(M_PI/(2*R)*(-R+n*h)))-1/(2*h*4*(-loesungsvektor[n-1]/h+1/R))*(h*pot((-R+(n+1)*h)/R, 4)-1/(5*(-loesungsvektor[n-1]/h+1/R))*(pot((-R+(n+1)*h)/R, 5)-pot(loesungsvektor[n-1]+(-R+n*h)/R, 5)))-1/h*(1/tau-(double) 1/2)*(pot(-R+(n+1)*h, 3)/(3*R)-((-R+n*h)*pot(-R+(n+1)*h, 2))/(2*R)+pot(-R+n*h, 3)/(6*R)))+1/tau*(-cos(M_PI/(2*R)*(-R+(n+1)*h))*2*R/M_PI+cos(M_PI/(2*R)*(-R+n*h))*2*R/M_PI)-1/(8*(-loesungsvektor[n-1]/h+1/R))*(pot((-R+(n+1)*h)/R, 4)-pot(loesungsvektor[n-1]+(-R+n*h)/R, 4))-(1/tau-(double) 1/2)*(pot(-R+(n+1)*h, 2)/(2*R)-pot(-R+n*h, 2)/(2*R))));

        //Cholesky-Verfahren
        LR(n, matrixDH);
        vwsubs(n, matrixDH, loesungsvektor, vektornegH);
        rwsubs(n, matrixDH, loesungsvektor);

        for(int i=0;i<n;i++)
        {
            loesungsvektor[i] = loesungsvektor[i] + kopieloesungsvektor[i];
        }
        k++;
    }

    //l(x) wird noch zur Loesung addiert
    for(int i=0;i<n;i++)
    {
        loesungsvektor[i] = loesungsvektor[i]+ (-R+(i+1)*h)/R;
    }
    ausgabe_Vektor(loesungsvektor, n);

    //Loesungsvektor wird in .txt Datei gespeichert
    FILE *file = fopen("loesungvarverfahren3.txt", "w");
    fprintf(file, "%5.8f ", -R);
    fprintf(file, "-1\n");
    for(int i=0;i<n;i++)
    {
        fprintf(file, "%5.8f ", -R+h*(i+1));
        fprintf(file, "%5.8f\n", loesungsvektor[i]);
    }
    fprintf(file, "%5.8f ", R);
    fprintf(file, "1\n");



}

//LR-Zerlegung
int LR(int n, double** A)
{
    if(A==NULL || n<1 || A[0][1]==0) return 1;
    for(int i=0;i<n-1;i++)
    {
        A[i][0] = A[i+1][0]/A[i][1];
        A[i+1][1] = A[i+1][1]-A[i][0]*A[i][2];
    }
	return 0;
}

//Vorwaertssubstitution
int vwsubs(int n, double** B, double* b, double* d)
{
	if(B == NULL || b == NULL) return 1;
	b[0] = d[0];
	for(int i=1;i<n;i++)
    {
        b[i] = d[i]-B[i-1][0]*b[i-1];
    }
	return 0;
}

//Rueckwaertssubstitution
int rwsubs(int n, double** B, double* b)
{
	if(B == NULL || b == NULL) return 1;
	b[n-1] = b[n-1]/B[n-1][1];
	for(int i=n-2;i>=0;i--)
	{
		b[i] = (b[i]-B[i][2]*b[i+1])/B[i][1];
	}
	return 0;
}

//Funktion zur Ausgabe von Vektoren
void ausgabe_Vektor(double* b, int n)
{
    for(int i=0;i<n;i++)
    {
        printf(" %5.16f",b[i]);
        printf("\n");
    }
    printf("\n");
}

//Funktion zur Ausgabe von Matrizen
void ausgabe_Matrix(double** matrix, int n, int m)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			printf(" %5.2f",matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

//Funktion zur Potenzrechnung
double pot(double d, int n)
{
    double ausgabe;
    ausgabe = 1;
    for(int i=0;i<n;i++)
    {
        ausgabe = ausgabe*d;
    }
    return ausgabe;
}



