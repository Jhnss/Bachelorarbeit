#include <stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

int LR(int n, double** A);
void ausgabe_Vektor(double* b, int n);
void ausgabe_Matrix(double** matrix, int n, int m);
int vwsubs(int n, double** B, double* b, double* d);
int rwsubs(int n, double** B, double* b);

int main(void)
{

    double** matrixDH;
	double* vektornegH;
	double* loesungsvektor;
	double* kopieloesungsvektor;
	int wiederholungen;//Anzahl der Widerholungen beim Newton-Verfahren
	int k;
	int n;// Anzahl der Zeitintervalle - 1
	double R;//Grenzen des Definitionsbereichs
	double h;//Schrittweite
	double tau;//Schrittweite der Diskretisierung des partiellen Diffgl
	wiederholungen = 10;
	n = 100;
	R = 5;
	h = (R+R)/(n+1);
	tau = 1.7;
    k = 0;

    //Anfangswert für Lösungsvektor wird erstellt
    loesungsvektor = (double*)calloc(n, sizeof(double));
    for(int i=0;i<n;i++)
	{
		loesungsvektor[i] = sin(M_PI/(2*R)*(-R+(i+1)*h)) ;
	}

    matrixDH = (double**)calloc(n, sizeof(double*));
    for(int i=0;i<n;i++)
    {
        matrixDH[i] = (double*)calloc(3, sizeof(double));
    }
    kopieloesungsvektor = (double*)calloc(n, sizeof(double));
    vektornegH = (double*)calloc(n, sizeof(double));

    //Ab hier beginnt das Newton Verfahren
	while(k<wiederholungen)
    {

        for(int i=0;i<n;i++)
            {
            kopieloesungsvektor[i] = loesungsvektor[i];
        }
        //DH wird erstellt
        for(int i=0;i<n;i++)
        {
            matrixDH[i][0] = -1;
            matrixDH[i][1] = 2+(1/tau-0.5)*h*h+h*h*3/2*loesungsvektor[i]*loesungsvektor[i];
            matrixDH[i][2] = -1;
        }
        matrixDH[0][0] = 0;
        matrixDH[n-1][2] = 0;
        //-H wird erstellt

        vektornegH[0] = -2*loesungsvektor[0]+loesungsvektor[1]-1+h*h*(((double) 1/2 - 1/tau)*loesungsvektor[0]+1/tau*sin(M_PI/(2*R)*(-R+h))- (double) 1/2*loesungsvektor[0]*loesungsvektor[0]*loesungsvektor[0]);
        for(int i=1;i<n-1;i++)
        {
        vektornegH[i] = -2*loesungsvektor[i]+loesungsvektor[i-1]+loesungsvektor[i+1]+h*h*(((double) 1/2 - 1/tau)*loesungsvektor[i]+1/tau*sin(M_PI/(2*R)*(-R+(i+1)*h))- (double) 1/2*loesungsvektor[i]*loesungsvektor[i]*loesungsvektor[i]);
        }
        vektornegH[n-1] = -2*loesungsvektor[n-1]+loesungsvektor[n-2]+1+h*h*(((double) 1/2 - 1/tau)*loesungsvektor[n-1]+1/tau*sin(M_PI/(2*R)*(-R+n*h))- (double) 1/2*loesungsvektor[n-1]*loesungsvektor[n-1]*loesungsvektor[n-1]);

        LR(n, matrixDH);
        vwsubs(n, matrixDH, loesungsvektor, vektornegH);
        rwsubs(n, matrixDH, loesungsvektor);
        for(int i=0;i<n;i++)
        {
            loesungsvektor[i] = loesungsvektor[i] + kopieloesungsvektor[i];
        }
        k++;
    }
    ausgabe_Vektor(loesungsvektor, n);

    //Lösungsvektor wird in .txt Datei gespeichert
    FILE *file = fopen("loesungdiffverfahren.txt", "w");
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
void ausgabe_Vektor(double* b, int n)
{
    for(int i=0;i<n;i++)
    {
        printf(" %5.8f",b[i]);
        printf("\n");
    }
    printf("\n");
}
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


