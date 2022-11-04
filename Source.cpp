
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include<cmath>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace cv;
using namespace std;


//Funcion para generar e kernell y normalizarlo
float *makeKernell(int n, float sigma) {
	int i, j;
	float *ker = (float*) malloc((n * n) * sizeof(float));
	float con =  1 / (2*M_PI*sigma*sigma) ;
	float s = (float)(n - 1) / 2;
	float x, y,sum=0.0;
	if (n % 2 == 0 || n<0) {
		printf("invlaid kernell size\n");
		return NULL;
	}


	else {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				x = j - s;
				y = i - s;
				*(ker + i*n + j) = con * exp(-(((x * x) + (y * y)) )/ (2 * sigma * sigma));
				sum += *(ker + i*n + j);

				printf("%f ", *(ker + i*n + j));
			}
			printf("\n");
		}
		
		printf("\n");
		float div = sum;
		float count = 0;
		
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				
				*(ker + i*n+ j) /= sum ;
				printf("%f ", *(ker + i * n + j));
				count += *(ker + i * n + j);
			}
			printf("\n");
		}
		printf("count: %f", count);


	}
	return ker;
}

// tamaño de la expansion de la matriz
int get_expand(int n) {
	if (n % 2 == 0) {
		return NULL;
	}
	return (n - 1) / 2;
}

int main() {
	//obtener los datos iniciales k y sigma
	int k_size;
	printf("tamaño kernell: ");
	scanf("%d", &k_size);
	float sigma = 0;
	printf("Sigma: ");
	scanf("%f", &sigma);
	float* kernell = makeKernell(k_size, sigma);
	int expan = get_expand(k_size);
	printf("Expansion size: %d\n", expan);

	//Leer imagen
	char NombreImagen[] = "lena.png";
	Mat image;
	image = imread(NombreImagen);
	int rows = image.rows;
	int columns = image.cols;
	Mat imageg(rows, columns, CV_8UC1);
	int i, j;
	double R, G, B, gray;


	if (!image.data)
	{
		cout << "\nError al cargar la imagen: " << NombreImagen << endl;
		exit(1);
	}
	//convertir la imagen RGB a grises
	for (i = 0; i < image.rows; i++) {
		for (j = 0; j < image.cols; j++) {

			R = image.at<Vec3b>(Point(j, i)).val[0];
			G = image.at<Vec3b>(Point(j, i)).val[1];
			B = image.at<Vec3b>(Point(j, i)).val[2];

			gray = R * 0.299 + G * 0.587 + B * 0.114;
			imageg.at<uchar>(Point(j, i)) = uchar(gray);
		}
	}


	//hacer la expansion a la matriz de grises para procesarla
	copyMakeBorder(imageg, imageg, expan, expan, expan, expan, BORDER_CONSTANT, 0);

	//inicializar el filtro de Gauss
	Mat gauss(image.rows, image.cols, CV_8UC1);
	int x = 0, y = 0;
	float sum = 0, pixel = 0;
	int k = k_size;



	//recorrer la matriz y asignar el valor del filtro de gauss a cada pixel
	for (i = expan; i < imageg.rows - expan; i++) {
		for (j = expan; j < imageg.cols - expan; j++) {


			for (y = 0; y < k; y++) {
				for (x = 0; x < k; x++) {

					sum += *(kernell + y * k + x) * (float)imageg.at<uchar>(Point(j + (x - expan), i + (y - expan)));

				}
			}

			gauss.at<uchar>(Point(j - expan, i - expan)) = uchar(floor(sum));
			sum = 0;

		}
	}


	
	


	//mostar imagenes obtenidas
	namedWindow("Original", WINDOW_AUTOSIZE);
	imshow("Original", image);

	namedWindow("Gris borde", WINDOW_AUTOSIZE);
	imshow("Gris borde", imageg);

	namedWindow("imagen suavizada", WINDOW_AUTOSIZE);
	imshow("imagen suavizada", gauss);

	//mostar dimensiones de cada imagen

	printf("\nImagen original");
	printf("\nrows: %d ", image.rows);
	printf("\ncolumns: %d \n", image.cols);

	printf("\nImagen en gris con borde");
	printf("\nrows: %d ", imageg.rows);
	printf("\ncolumns: %d \n", imageg.cols);

	printf("\nImagen suavizada");
	printf("\nrows: %d ", gauss.rows);
	printf("\ncolumns: %d \n", gauss.cols);

	//ecualizar la imagen con el filtro de gauss
	equalizeHist(gauss, gauss);
	namedWindow("imagen equalizada", WINDOW_AUTOSIZE);
	imshow("imagen equalizada", gauss);


	printf("\nIimagen equlizada");
	printf("\nrows: %d ", gauss.rows);
	printf("\ncolumns: %d \n", gauss.cols);


	//inicializar las variables requeridas para obtener |G| y deteccion de bordes
	Mat temp(rows,columns, CV_8UC1);
	Mat Fx(rows, columns, CV_8UC1);
	Mat Fy(rows, columns, CV_8UC1);
	Mat Gval(rows, columns, CV_8UC1);
	Mat Gori(rows, columns, CV_8UC1);
	Mat line(rows, columns, CV_8UC1);

	//expandir las matrizes sobre las que se van a trabajar
	copyMakeBorder(gauss, temp, 1, 1, 1, 1, BORDER_CONSTANT, 0);
	copyMakeBorder(gauss, Gori, 1, 1, 1, 1, BORDER_CONSTANT, 0);

	//kernells que se van a aplicar
	int hx[3][3] = { -1, 0, 1 , -2, 0, 2, -1, 0, 1 };
	int hy[3][3] = { 1, 2, 1 , 0, 0, 0, -1, -2, -1 };

	
	float sumx = 0, sumy = 0;
	int fx,fy;
 int theta;
 int ax, ay;
 float max=0;
 int ybot[2], ytop[2], xest;
 float tmp;

 //recorrer la matriz equalizada y con el filtro de gauss y hacer las operaciones correspondientes
	for (i = 1; i < temp.rows-1; i++) {
		for (j = 1; j < temp.cols - 1; j++) {


			//aplicar los kernells hx y hy a sus respectivas matrizes
			for (y = 0; y < 3; y++) {
				for (x = 0; x < 3; x++) {
					sumx += hx[y][x] *(float) temp.at<uchar>(Point(j + (x - 1), i + (y - 1)));
					sumy += hy[y][x] *(float) temp.at<uchar>(Point(j + (x - 1), i + (y - 1)));
				}
			}

			//si el valor del pixel obtenido es menor a cero se descarta el pixel
			if (sumx < 0) {
				sumx = 0;
			}
			if (sumy < 0) {
				sumy = 0;
			}

			//llenamos las matrices Fx y Fy
			Fx.at<uchar>(Point(j-1, i-1)) = floor(abs(sumx));
			Fy.at<uchar>(Point(j-1, i-1)) = floor(abs(sumy));

			fx = Fx.at<uchar>(Point(j - 1, i - 1));
			fy = Fy.at<uchar>(Point(j - 1, i - 1));
			 
			//obtenemos la orientacion de G si un pixel es 0 se aproxima el valor lo mas posible
			if (fx == 0) {
				Gori.at<uchar>(Point(j - 1, i - 1)) = round(atan(fy / 0.000000001) * (180 / M_1_PI));
			}
			else {
				Gori.at<uchar>(Point(j - 1, i - 1)) = round(atan(fy / fx) * (180 / M_1_PI));
			}
			//obtenemos |G| y la copiamos para aplicar la etccion de lineas
			Gval.at<uchar>(Point(j - 1, i - 1))=  floor(sqrt(fx*fx+fy*fy));
			line.at<uchar>(Point(j - 1, i - 1)) = Gval.at<uchar>(Point(j - 1, i - 1));
			theta = Gori.at<uchar>(Point(j - 1, i - 1));

			//Supresion de no maximos sin interpolacion, se descarto porque el otro metodo daba
			//mejores resultados
			/*line.at<uchar>(Point(j - 1, i - 1)) = Gori.at<uchar>(Point(j - 1, i - 1));
			if ((theta >= 337) || (theta < 22)) { ax = 1; ay = 0; }
			else if ((theta >= 22.5) && (theta < 67.5)) { ax = 1 ; ay = 1; }
			else if ((theta >= 67.5) && (theta < 112.5)) { ax = 0; ay = 1; }
			else if ((theta >= 112.5) && (theta < 157.5)) { ax = -1; ay = 1; }
			else if ((theta >= 157.5) && (theta < 202.5)) { ax = -1; ay = 0; }
			else if ((theta >= 202.5) && (theta < 247.5)) { ax = -1; ay = -1; }
			else if ((theta >= 247.5) && (theta < 292.5)) { ax = 0; ay = -1; }
			else if ((theta >= 292.5) && (theta < 337.5)) { ax = 1; ay = -1; }
			if (theta < Gori.at<uchar>(Point(j + (-ax), i + (-ay))) || theta < Gori.at<uchar>(Point(j + ax, i + ay))) {
				Gori.at<uchar>(Point(j - 1, i - 1)) = 0;}*/

			//supresion de no maximos con interpolacion
			if (Gval.at<uchar>(Point(j - 1, i - 1)) == 0) { Gval.at<uchar>(Point(j - 1, i - 1)) = 10; }
			if ((theta >= 0 && theta <= 45) || (theta < 225 && theta >= 270)) {
				ybot[0] = Gori.at<uchar>(Point(j + 1, i));
				ybot[1] = Gori.at<uchar>(Point(j + 1, i + 1));

				ytop[0] = Gori.at<uchar>(Point(j - 1, i));
				ytop[1] = Gori.at<uchar>(Point(j - 1, i - 1));

				xest = abs(fy / Gval.at<uchar>(Point(j - 1, i - 1)));
				if(theta< ((ybot[1]-ybot[0])*xest+ybot[0]) && theta < ((ytop[1] - ytop[0]) * xest + ytop[0]))
					line.at<uchar>(Point(j - 1, i - 1)) = 0;
			}
			else if ((theta > 45 && theta <= 90) || (theta < 270 && theta >= 225)) {
				ybot[0] = Gori.at<uchar>(Point(j, i+1));
				ybot[1] = Gori.at<uchar>(Point(j + 1, i + 1));

				ytop[0] = Gori.at<uchar>(Point(j , i-1));
				ytop[1] = Gori.at<uchar>(Point(j - 1, i - 1));

				xest = abs(fy / Gval.at<uchar>(Point(j - 1, i - 1)));
				if (theta < ((ybot[1] - ybot[0]) * xest + ybot[0]) && theta < ((ytop[1] - ytop[0]) * xest + ytop[0]))
					line.at<uchar>(Point(j - 1, i - 1)) = 0;
			}
			else if ((theta > 90 && theta <= 135) || (theta < 315 && theta >= 270)) {
				ybot[0] = Gori.at<uchar>(Point(j , i+1));
				ybot[1] = Gori.at<uchar>(Point(j - 1, i + 1));

				ytop[0] = Gori.at<uchar>(Point(j, i-1));
				ytop[1] = Gori.at<uchar>(Point(j +1 , i - 1));

				xest = abs(fy / Gval.at<uchar>(Point(j - 1, i - 1)));
				if (theta < ((ybot[1] - ybot[0]) * xest + ybot[0]) && theta < ((ytop[1] - ytop[0]) * xest + ytop[0]))
					line.at<uchar>(Point(j - 1, i - 1)) = 0;
			}
			else if ((theta > 135 && theta <= 180) || (theta < 360 && theta >= 315)) {
				ybot[0] = Gori.at<uchar>(Point(j - 1, i));
				ybot[1] = Gori.at<uchar>(Point(j - 1, i + 1));

				ytop[0] = Gori.at<uchar>(Point(j +1, i));
				ytop[1] = Gori.at<uchar>(Point(j + 1, i - 1));

				xest = abs(fy / Gval.at<uchar>(Point(j - 1, i - 1)));
				if (theta < ((ybot[1] - ybot[0]) * xest + ybot[0]) && theta < ((ytop[1] - ytop[0]) * xest + ytop[0]))
					line.at<uchar>(Point(j - 1, i - 1)) = 0;
			}

			
				
		
			
			//buscamos el maximo valor
			if (max < line.at<uchar>(Point(j - 1, i - 1))) {
				max = line.at<uchar>(Point(j - 1, i - 1));
			}

			sumx = 0;
			sumy = 0;

		}
	}

	//aplicamos los limites
for (i = 1; i < temp.rows - 1; i++) {
		for (j = 1; j < temp.cols - 1; j++) {

			//limite usando el metodo sin interpolacion
			/*tmp = Gori.at<uchar>(Point(j - 1, i - 1));
			if (tmp > 90 && tmp<=180) { tmp = (180 - tmp); 
			}
			else if (tmp > 180 && 270 >= tmp) { tmp = (270 - tmp);
			}
			else if ( tmp > 270 && 360 >= tmp) { tmp = (360 - tmp);
			}




			if(tmp<60)
				Gori.at<uchar>(Point(j - 1, i - 1))=0;
			if (tmp > 75) {
				Gori.at<uchar>(Point(j - 1, i - 1)) = 255;
			}*/

			//metodo con interpolacion borde descartable < 15% borde fuerte > 80%
			tmp = (float)line.at<uchar>(Point(j - 1, i - 1))/max;
			if (tmp < .15)
				line.at<uchar>(Point(j - 1, i - 1)) = 0;
			if (tmp > .80) {
				line.at<uchar>(Point(j - 1, i - 1)) = 255;
			}
				

		}
	}


	//mostar resultados
	namedWindow("|G|", WINDOW_AUTOSIZE);
	imshow("|G|", Gval);


	namedWindow("Deteccion de contornos", WINDOW_AUTOSIZE);
	imshow("Deteccion de contornos", line);

	printf("\n|G|");
	printf("\nrows: %d ", Gval.rows);
	printf("\ncolumns: %d \n", Gval.cols);

	printf("\nDeteccion de bordes");
	printf("\nrows: %d ", line.rows);
	printf("\ncolumns: %d \n", line.cols);

	waitKey(0);
	return 1;

}