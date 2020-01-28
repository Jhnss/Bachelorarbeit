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
	double* vektorNegH;
	double* vektorUk;
	double* vektorLoesung;
	double* vektorLoesungKopie;
    int methode;//Zur Auswahl der Methode
    int wiederholungen;//Anzahl der Wiederholungen beim Newton-Verfahren
    int zeitschritte;//Anzahl der Zeitschritte
	int n;// Anzahl der Intervalle, in die [-R,R] geteilt werden soll - 1
	double R;//Grenzen des Definitionsbereichs
	double h;//Länge der Intervalle, in die [-R,R] unterteilt wird
	double tau;//Schrittweite der Diskretisierung des partiellen Diffgl

	//Einstellung der Parameter
	printf("Waehlen Sie zwischen Variationsmethode(0) oder Differenzenverfahren(1)\n");
	scanf("%d", &methode);
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
    vektorUk = (double*)calloc(n, sizeof(double));
    matrixDH = (double**)calloc(n, sizeof(double*));
    for(int i=0;i<n;i++)
    {
        matrixDH[i] = (double*)calloc(3, sizeof(double));
    }
    vektorLoesungKopie = (double*)calloc(n, sizeof(double));
    vektorNegH = (double*)calloc(n, sizeof(double));

    if(methode==1) differenzenverfahren(matrixDH, vektorNegH, vektorLoesung, vektorLoesungKopie, vektorUk, wiederholungen, zeitschritte, n, R, h, tau);
    else variationsmethode(matrixDH, vektorNegH, vektorLoesung, vektorLoesungKopie, vektorUk, wiederholungen, zeitschritte, n, R, h, tau);
}
int differenzenverfahren(double** matrixDH, double* vektorNegH, double* vektorLoesung, double* vektorLoesungKopie, double* vektorUk, int wiederholungen, int zeitschritte, int n, double R, double h, double tau)
{
    //Anfangswert fuer das Newton_Verfahren wird festgelegt
    for(int i=0;i<n;i++)
	{
		vektorLoesung[i] = 1/5*(-R+(i+1)*h);
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

int variationsmethode(double** matrixDH, double* vektorNegH, double* vektorLoesung, double* vektorLoesungKopie, double* vektorUk, int wiederholungen, int zeitschritte, int n, double R, double h, double tau)
{
    double matrixA_diagnonaleintrag;
    double matrixA_nichtDiagonaleintrag;

    //Anfangswert fuer das Newton_Verfahren wird festgelegt
    for(int i=0;i<n;i++)
	{
		vektorLoesung[i] = 1/5*(-R+(i+1)*h)- (-R+(i+1)*h)/R;
	}

    //Berechnung der Matrix A
	matrixA_diagnonaleintrag = 2/h+(1/tau-0.5)*2*h/3;
	matrixA_nichtDiagonaleintrag = -1/h+(1/tau-0.5)*1/6*h;

    //Anfang der Berechnungen der Loesung
    for(int k=0;k<zeitschritte;k++)
    {
        //Vektor u_k wird erstellt
        for(int i=0;i<n;i++)
        {
            vektorUk[i] = vektorLoesung[i]+ (-R+(i+1)*h)/R;
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
                if(i>0 && i<(n-1)) matrixDH[i][1] = matrixA_diagnonaleintrag+3/2*(1/(3*(-vektorLoesung[i-1]/h+vektorLoesung[i]/h+1/R))*pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 3)-2/(3*h*(-vektorLoesung[i-1]/h+vektorLoesung[i]/h+1/R))*(1/(4*(-vektorLoesung[i-1]/h+vektorLoesung[i]/h+1/R))*pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 4)-1/(20*h*pot(-vektorLoesung[i-1]/h+vektorLoesung[i]/h+1/R, 2))*(pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 5)-pot(vektorLoesung[i-1]+(-R+i*h)/R, 5)))-1/(3*(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R))*pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 3)+2/(3*h*(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R))*(-1/(4*(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R))*pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 4)+1/(20*h*pot(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R, 2))*(pot(vektorLoesung[i+1]+(-R+(i+2)*h)/R, 5)-pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 5))));
                if(i<(n-1)) matrixDH[i][2] = matrixA_nichtDiagonaleintrag-(1/(2*(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R))*(-1/(h*4*(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R))*(pot(vektorLoesung[i+1]+(-R+(i+2)*h)/R, 4)+pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 4))+1/(10*h*h*pot(-vektorLoesung[i]/h+vektorLoesung[i+1]/h+1/R, 2))*(pot(vektorLoesung[i+1]+(-R+(i+2)*h)/R, 5)-pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 5))));
                if(i>0) matrixDH[i][0] = matrixDH[i-1][2];
            }
            matrixDH[0][1] = matrixA_diagnonaleintrag+3/2*(1/(3*(vektorLoesung[0]/h+1/R))*pot(vektorLoesung[0]+(-R+h)/R, 3)-2/(3*h*(vektorLoesung[0]/h+1/R))*(1/(4*(vektorLoesung[0]/h+1/R))*pot(vektorLoesung[0]+(-R+h)/R, 4)-1/(20*h*pot(vektorLoesung[0]/h+1/R, 2))*(pot(vektorLoesung[0]+(-R+h)/R, 5)+1))-1/(3*(-vektorLoesung[0]/h+vektorLoesung[1]/h+1/R))*pot(vektorLoesung[0]+(-R+h)/R, 3)+2/(3*h*(-vektorLoesung[0]/h+vektorLoesung[1]/h+1/R))*(-1/(4*(-vektorLoesung[0]/h+vektorLoesung[1]/h+1/R))*pot(vektorLoesung[0]+(-R+h)/R, 4)+1/(20*h*pot(-vektorLoesung[0]/h+vektorLoesung[1]/h+1/R, 2))*(pot(vektorLoesung[1]+(-R+2*h)/R, 5)-pot(vektorLoesung[0]+(-R+h)/R, 5))));
            matrixDH[n-1][1] = matrixA_diagnonaleintrag+3/2*(1/(3*(-vektorLoesung[n-2]/h+vektorLoesung[n-1]/h+1/R))*pot(vektorLoesung[n-1]+(-R+n*h)/R, 3)-2/(3*h*(-vektorLoesung[n-2]/h+vektorLoesung[n-1]/h+1/R))*(1/(4*(-vektorLoesung[n-2]/h+vektorLoesung[n-1]/h+1/R))*pot(vektorLoesung[n-1]+(-R+n*h)/R, 4)-1/(20*h*pot(-vektorLoesung[n-2]/h+vektorLoesung[n-1]/h+1/R, 2))*(pot(vektorLoesung[n-1]+(-R+n*h)/R, 5)-pot(vektorLoesung[n-2]+(-R+(n-1)*h)/R, 5)))-1/(3*(-vektorLoesung[n-1]/h+1/R))*pot(vektorLoesung[n-1]+(-R+n*h)/R, 3)+2/(3*h*(-vektorLoesung[n-1]/h+1/R))*(-1/(4*(-vektorLoesung[n-1]/h+1/R))*pot(vektorLoesung[n-1]+(-R+n*h)/R, 4)+1/(20*h*pot(-vektorLoesung[n-1]/h+1/R, 2))*(1-pot(vektorLoesung[n-1]+(-R+n*h)/R, 5))));

            //-H wird erstellt
            vektorNegH[0] = -(vektorLoesung[0]*matrixA_diagnonaleintrag+matrixA_nichtDiagonaleintrag*vektorLoesung[1]-(1/(tau*2*h*(vektorUk[0]/h+1/h))*(h*pot(vektorUk[0], 2)-1/(3*(vektorUk[0]/h+1/h))*(pot(vektorUk[0], 3)-pot(-1, 3)))-1/(2*h*4*(vektorLoesung[0]/h+1/R))*(h*pot(vektorLoesung[0]+(-R+h)/R, 4)-1/(5*(vektorLoesung[0]/h+1/R))*(pot(vektorLoesung[0]+(-R+h)/R, 5)-pot(-R/R, 5)))-1/h*(1/tau-0.5)*(pot(-R+h, 3)/(3*R)-((-R)*pot(-R+h, 2))/(2*R)+pot(-R, 3)/(6*R)))-(-1*(1/(tau*2*h*(vektorUk[1]/h-vektorUk[0]/h))*(h*pot(vektorUk[1], 2)-1/(3*(vektorUk[1]/h-vektorUk[0]/h))*(pot(vektorUk[1], 3)-pot(vektorUk[0], 3)))-1/(2*h*4*(vektorLoesung[1]/h-vektorLoesung[0]/h+1/R))*(h*pot(vektorLoesung[1]+(-R+2*h)/R, 4)-1/(5*(vektorLoesung[1]/h-vektorLoesung[0]/h+1/R))*(pot(vektorLoesung[1]+(-R+2*h)/R, 5)-pot(vektorLoesung[0]+(-R+h)/R, 5)))-1/h*(1/tau-0.5)*(pot(-R+2*h, 3)/(3*R)-((-R+h)*pot(-R+2*h, 2))/(2*R)+pot(-R+h, 3)/(6*R)))+1/(2*tau*(vektorUk[1]/h-vektorUk[0]/h))*(pot(vektorUk[1], 2)-pot(vektorUk[0], 2))-1/(8*(vektorLoesung[1]/h-vektorLoesung[0]/h+1/R))*(pot(vektorLoesung[1]+(-R+2*h)/R, 4)-pot(vektorLoesung[0]+(-R+h)/R, 4)) -(1/tau-0.5)*(pot(-R+2*h, 2)/(2*R)-pot(-R+h, 2)/(2*R))));
            for(int i=1;i<n-1;i++)
            {
                vektorNegH[i] = -(vektorLoesung[i]*matrixA_diagnonaleintrag+matrixA_nichtDiagonaleintrag*(vektorLoesung[i-1]+vektorLoesung[i+1])-(1/(tau*2*h*(vektorUk[i]/h-vektorUk[i-1]/h))*(h*pot(vektorUk[i], 2)-1/(3*(vektorUk[i]/h-vektorUk[i-1]/h))*(pot(vektorUk[i], 3)-pot(vektorUk[i-1], 3)))-1/(2*h*4*(vektorLoesung[i]/h-vektorLoesung[i-1]/h+1/R))*(h*pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 4)-1/(5*(vektorLoesung[i]/h-vektorLoesung[i-1]/h+1/R))*(pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 5)-pot(vektorLoesung[i-1]+(-R+i*h)/R, 5)))-1/h*(1/tau-0.5)*(pot(-R+(i+1)*h, 3)/(3*R)-((-R+i*h)*pot(-R+(i+1)*h, 2))/(2*R)+pot(-R+i*h, 3)/(6*R)))-(-1*(1/(tau*2*h*(vektorUk[i+1]/h-vektorUk[i]/h))*(h*pot(vektorUk[i+1], 2)-1/(3*(vektorUk[i+1]/h-vektorUk[i]/h))*(pot(vektorUk[i+1], 3)-pot(vektorUk[i], 3)))-1/(2*h*4*(vektorLoesung[i+1]/h-vektorLoesung[i]/h+1/R))*(h*pot(vektorLoesung[i+1]+(-R+(i+2)*h)/R, 4)-1/(5*(vektorLoesung[i+1]/h-vektorLoesung[i]/h+1/R))*(pot(vektorLoesung[i+1]+(-R+(i+2)*h)/R, 5)-pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 5)))-1/h*(1/tau-0.5)*(pot(-R+(i+2)*h, 3)/(3*R)-((-R+(i+1)*h)*pot(-R+(i+2)*h, 2))/(2*R)+pot(-R+(i+1)*h, 3)/(6*R)))+1/(2*tau*(vektorUk[i+1]/h-vektorUk[i]/h))*(pot(vektorUk[i+1], 2)-pot(vektorUk[i], 2))-1/(8*(vektorLoesung[i+1]/h-vektorLoesung[i]/h+1/R))*(pot(vektorLoesung[i+1]+(-R+(i+2)*h)/R, 4)-pot(vektorLoesung[i]+(-R+(i+1)*h)/R, 4))-(1/tau-0.5)*(pot(-R+(i+2)*h, 2)/(2*R)-pot(-R+(i+1)*h, 2)/(2*R))));
            }
            vektorNegH[n-1] = -(vektorLoesung[n-1]*matrixA_diagnonaleintrag+matrixA_nichtDiagonaleintrag*vektorLoesung[n-2]-(1/(tau*2*h*(vektorUk[n-1]/h-vektorUk[n-2]/h))*(h*pot(vektorUk[n-1], 2)-1/(3*(vektorUk[n-1]/h-vektorUk[n-2]/h))*(pot(vektorUk[n-1], 3)-pot(vektorUk[n-2], 3)))-1/(2*h*4*(vektorLoesung[n-1]/h-vektorLoesung[n-2]/h+1/R))*(h*pot(vektorLoesung[n-1]+(-R+n*h)/R, 4)-1/(5*(vektorLoesung[n-1]/h-vektorLoesung[n-2]/h+1/R))*(pot(vektorLoesung[n-1]+(-R+n*h)/R, 5)-pot(vektorLoesung[n-2]+(-R+(n-1)*h)/R, 5)))-1/h*(1/tau-0.5)*(pot(-R+n*h, 3)/(3*R)-((-R+(n-1)*h)*pot(-R+n*h, 2))/(2*R)+pot(-R+(n-1)*h, 3)/(6*R)))-(-1*(1/(tau*2*h*(1/h-vektorUk[n-1]/h))*(h*pot(1, 2)-1/(3*(1/h-vektorUk[n-1]/h))*(pot(1, 3)-pot(vektorUk[n-1], 3)))-1/(2*h*4*(-vektorLoesung[n-1]/h+1/R))*(h*pot((-R+(n+1)*h)/R, 4)-1/(5*(-vektorLoesung[n-1]/h+1/R))*(pot((-R+(n+1)*h)/R, 5)-pot(vektorLoesung[n-1]+(-R+n*h)/R, 5)))-1/h*(1/tau-0.5)*(pot(-R+(n+1)*h, 3)/(3*R)-((-R+n*h)*pot(-R+(n+1)*h, 2))/(2*R)+pot(-R+n*h, 3)/(6*R)))+1/(2*tau*(1/h-vektorUk[n-1]/h))*(pot(1, 2)-pot(vektorUk[n-1], 2))-1/(8*(-vektorLoesung[n-1]/h+1/R))*(pot((-R+(n+1)*h)/R, 4)-pot(vektorLoesung[n-1]+(-R+n*h)/R, 4))-(1/tau-0.5)*(pot(-R+(n+1)*h, 2)/(2*R)-pot(-R+n*h, 2)/(2*R))));

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

    //l(x) wird noch zur Loesung addiert
    for(int i=0;i<n;i++)
    {
        vektorLoesung[i] = vektorLoesung[i]+ (-R+(i+1)*h)/R;
    }

    printf("\nDer Loesungsvektor ist:\n");
    ausgabe_Vektor(vektorLoesung, n);

    //Loesungsvektor wird in .txt Datei gespeichert
    FILE *file = fopen("loesungVariationsmethode.txt", "w");
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



