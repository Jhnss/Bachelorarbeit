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
	double* vektorNegH;
	double* vektorLoesung;
	double* vektorLoesungKopie;
	double* vektorUk;
	int wiederholungen;//Anzahl der Wiederholungen beim Newton-Verfahren
    int zeitschritte;//Anzahl der Zeitschritte
	int n;// Anzahl der Intervalle, in die [-R,R] geteilt werden soll - 1
	double R;//Grenzen des Definitionsbereichs
	double h;//Länge der Intervalle, in die [-R,R] unterteilt wird
	double tau;//Schrittweite der Diskretisierung des partiellen Diffgl

	//Einstellung der Parameter
	printf("Differenzenverfahren\n");
	printf("Bestimmen Sie den Parameter R:\n");
	scanf("%lf", &R);
	printf("Bestimmen Sie den Zeitschritt tau:\n");
	scanf("%lf", &tau);
	printf("Bestimmen Sie die Anzahl der Zeitschritte:\n");
	scanf("%d", &zeitschritte);
	printf("Bestimmen Sie die Anzahl der Intervalle, in die [-R,R] geteilt werden soll:\n");
	scanf("%d", &n);
    printf("Bestimmen Sie die Anzahl der Wiederholungen beim Newton Verfahren:\n");
	scanf("%d", &wiederholungen);
	n = n-1;
	h = (R+R)/(n+1);

    //Speicherplatz für Matrizen und Vektoren wird festgelegt
    vektorLoesung = (double*)calloc(n, sizeof(double));
    matrixDH = (double**)calloc(n, sizeof(double*));
    for(int i=0;i<n;i++)
    {
        matrixDH[i] = (double*)calloc(3, sizeof(double));
    }

    vektorLoesungKopie = (double*)calloc(n, sizeof(double));
    vektorNegH = (double*)calloc(n, sizeof(double));
    vektorUk = (double*)calloc(n, sizeof(double));

    //Anfangswert fuer das Newton_Verfahren wird festgelegt
    for(int i=0;i<n;i++)
	{
		vektorLoesung[i] = sin(M_PI/(2*R)*(-R+(i+1)*h));
	}

    //Anfang der Berechnungen der Loesung
    for(int k=0;k<zeitschritte;k++)
    {
        //Vektor u_k wird erstellt
        for(int i=0;i<n;i++)
        {
            vektorUk[i] = vektorLoesung[i];
        }

        //Beginn des Newton-Verfahrens
        for(int j=0;j<wiederholungen;j++)
        {
            //Kopie des Loesungsvektor wird erstellt
            for(int i=0;i<n;i++)
            {
                vektorLoesungKopie[i] = vektorLoesung[i];
            }

            //DH wird erstellt
            for(int i=0;i<n;i++)
            {
                matrixDH[i][0] = -1;
                matrixDH[i][1] = 2+(1/tau-0.5)*h*h+h*h*3/2*vektorLoesung[i]*vektorLoesung[i];
                matrixDH[i][2] = -1;
            }

            //-H wird erstellt
            vektorNegH[0] = -2*vektorLoesung[0]+vektorLoesung[1]-1+h*h*((0.5 - 1/tau)*vektorLoesung[0]+1/tau*vektorUk[0]-0.5*vektorLoesung[0]*vektorLoesung[0]*vektorLoesung[0]);
            for(int i=1;i<n-1;i++)
            {
                vektorNegH[i] = -2*vektorLoesung[i]+vektorLoesung[i-1]+vektorLoesung[i+1]+h*h*((0.5 - 1/tau)*vektorLoesung[i]+1/tau*vektorUk[i]-0.5*vektorLoesung[i]*vektorLoesung[i]*vektorLoesung[i]);
            }
            vektorNegH[n-1] = -2*vektorLoesung[n-1]+vektorLoesung[n-2]+1+h*h*((0.5 - 1/tau)*vektorLoesung[n-1]+1/tau*vektorUk[n-1]-0.5*vektorLoesung[n-1]*vektorLoesung[n-1]*vektorLoesung[n-1]);

            //Beginn des Cholesky-Verfahrens
            LR(n, matrixDH);
            vwsubs(n, matrixDH, vektorLoesung, vektorNegH);
            rwsubs(n, matrixDH, vektorLoesung);

            for(int i=0;i<n;i++)
            {
                vektorLoesung[i] = vektorLoesung[i] + vektorLoesungKopie[i];
            }
        }
    }

    printf("\nDer Loesungsvektor ist:\n");
    ausgabe_Vektor(vektorLoesung, n);

    //Loesungsvektor wird in .txt Datei gespeichert
    FILE *file = fopen("loesungDifferenzenverfahren.txt", "w");
    fprintf(file, "%5.8f ", -R);
    fprintf(file, "-1\n");
    for(int i=0;i<n;i++)
    {
        fprintf(file, "%5.8f ", -R+h*(i+1));
        fprintf(file, "%5.8f\n", vektorLoesung[i]);
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



