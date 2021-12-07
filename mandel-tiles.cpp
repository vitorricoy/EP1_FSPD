#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <queue>
#include <math.h>
#include <time.h>

#define MAXX 640
#define MAXY 480
#define MAXITER 32768

// params for each call to the fractal function
typedef struct {
	int left; int low;  // lower left corner in the screen
	int ires; int jres; // resolution in pixels of the area to compute
	double xmin; double ymin;   // lower left corner in domain (x,y)
	double xmax; double ymax;   // upper right corner in domain (x,y)
} fractal_param_t;

FILE* entrada; // descriptor for the list of tiles (cannot be stdin)
std::queue<fractal_param_t> tarefas;
const int numeroThreads = 5;
pthread_mutex_t acessoVariaveisEstatistica;
pthread_mutex_t acessoListaTarefas;
pthread_cond_t preencherTarefas;
pthread_cond_t esperarTarefas;
int totalTarefasExecutadas = 0;
int tarefasExecutadasThread[numeroThreads];
double tempoTotalGasto = 0.0;
double tempoTarefasThread[numeroThreads];
int contadorFilaVazia = 0;

void* threadTrabalhadora(void* param) {
	long threadId = (long) param;
	bool leuEOW = false;
	while(!leuEOW) {
		bool filaFicouVazia = false;
		pthread_mutex_lock(&acessoListaTarefas);
		while(tarefas.empty()) {
			filaFicouVazia = true;
			pthread_cond_wait(&esperarTarefas, &acessoListaTarefas);
		}
		fractal_param_t p =	tarefas.front();
		tarefas.pop();
		if(tarefas.size() == numeroThreads) {
			pthread_cond_signal(&preencherTarefas);
		}
		pthread_mutex_unlock(&acessoListaTarefas);
		clock_t inicio = clock();
		if(p.jres == 0 && p.ires == 0) {
			leuEOW = true;
		}
		if(!leuEOW) {
			double dx, dy;
			int i, j, k, color;
			double x, y, u, v, u2, v2;

			dx = ( p.xmax - p.xmin ) / p.ires;
			dy = ( p.ymax - p.ymin ) / p.jres;
			
			// scanning every point in that rectangular area.
			// Each point represents a Complex number (x + yi).
			// Iterate that complex number
			for(j = 0; j < p.jres; j++) {
				for(i = 0; i <= p.ires; i++) {
					x = i * dx + p.xmin; // c_real
					u = u2 = 0; // z_real
					y = j * dy + p.ymin; // c_imaginary
					v = v2 = 0; // z_imaginary

					// Calculate whether c(c_real + c_imaginary) belongs
					// to the Mandelbrot set or not and draw a pixel
					// at coordinates (i, j) accordingly
					// If you reach the Maximum number of iterations
					// and If the distance from the origin is
					// greater than 2 exit the loop
					for(k=0; (k < MAXITER) && ((u2+v2) < 4); ++k) {
						// Calculate Mandelbrot function
						// z = z*z + c where z is a complex number

						// imag = 2*z_real*z_imaginary + c_imaginary
						v = 2 * u * v + y;
						// real = z_real^2 - z_imaginary^2 + c_real
						u  = u2 - v2 + x;
						u2 = u * u;
						v2 = v * v;
					}
				}
			}
			clock_t fim = clock();
			double tempoGasto = ((double)fim - inicio) / CLOCKS_PER_SEC * 1000;
			pthread_mutex_lock(&acessoVariaveisEstatistica);
			totalTarefasExecutadas++;
			tarefasExecutadasThread[threadId]++;
			tempoTotalGasto += tempoGasto;
			tempoTarefasThread[threadId] += tempoGasto;
			if(filaFicouVazia) {
				contadorFilaVazia++;
			}
			pthread_mutex_unlock(&acessoVariaveisEstatistica);
		}
	}
	pthread_exit(NULL);
}

void* threadLeitora(void* param) {
	fractal_param_t p;
	int n;
	pthread_mutex_lock(&acessoListaTarefas);
	bool acabouArquivo = false;
	n = fscanf(entrada,"%d %d %d %d",&(p.left),&(p.low),&(p.ires),&(p.jres));
	if (n == EOF) acabouArquivo = true;
	while(!acabouArquivo) {
		if(n!=4) {
			perror("fscanf(left,low,ires,jres)");
			exit(-1);
		}
		n = fscanf(entrada,"%lf %lf %lf %lf",
			&(p.xmin),&(p.ymin),&(p.xmax),&(p.ymax));
		if(n!=4) {
			perror("scanf(xmin,ymin,xmax,ymax)");
			exit(-1);
		}
		
		tarefas.push(p);
		if(tarefas.size() >= 4*numeroThreads) {
			pthread_cond_broadcast(&esperarTarefas);
			pthread_cond_wait(&preencherTarefas, &acessoListaTarefas);
		}
		n = fscanf(entrada,"%d %d %d %d",&(p.left),&(p.low),&(p.ires),&(p.jres));
		if(n == EOF) acabouArquivo = true;
	}
	p.ires = 0;
	p.jres = 0;
	for(int i=0; i < numeroThreads; i++) {
		tarefas.push(p);
	}
	pthread_mutex_unlock(&acessoListaTarefas);
	pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
	int i,j,k;
	fractal_param_t p;

	if(argc!=2) {
		fprintf(stderr,"usage %s filename\n", argv[0] );
		exit(-1);
	}
	
	if((entrada=fopen(argv[1],"r"))==NULL) {
		perror("fdopen");
		exit(-1);
	}

	pthread_mutex_init(&acessoVariaveisEstatistica, NULL);
	pthread_mutex_init(&acessoListaTarefas, NULL);
	pthread_cond_init(&preencherTarefas, NULL);
	pthread_cond_init(&esperarTarefas, NULL);

	// criar threads
	pthread_t threadLeitura;
	pthread_create(&threadLeitura, NULL, threadLeitora, NULL);

	pthread_t threadsTrabalhadoras[numeroThreads];
	for(long i=0; i<numeroThreads; i++) {
		pthread_create(&threadsTrabalhadoras[i], NULL, threadTrabalhadora, (void*) i);
	}

	pthread_join(threadLeitura, NULL);

	for(int i=0; i<numeroThreads; i++) {
		pthread_join(threadsTrabalhadoras[i], NULL);
	}

	double mediaTarefas = (double) totalTarefasExecutadas/(double) numeroThreads;
	double desvioPadraoTarefas = 0;
	for(int i=0; i<numeroThreads; i++) {
		desvioPadraoTarefas += (tarefasExecutadasThread[i] - mediaTarefas)*(tarefasExecutadasThread[i] - mediaTarefas);
	}
	desvioPadraoTarefas /= (double) numeroThreads;
	desvioPadraoTarefas = sqrt(desvioPadraoTarefas);

	double mediaTempo = tempoTotalGasto / (double) numeroThreads;
	double desvioPadraoTempoGasto = 0;
	for(int i=0; i<numeroThreads; i++) {
		desvioPadraoTempoGasto += (tempoTarefasThread[i] - mediaTarefas)*(tempoTarefasThread[i] - mediaTarefas);
	}
	desvioPadraoTempoGasto /= (double) numeroThreads;
	desvioPadraoTempoGasto = sqrt(desvioPadraoTempoGasto);
	printf("Tarefas: total = %d;  média por trabalhador = %lf(%lf)\n", totalTarefasExecutadas, mediaTarefas, desvioPadraoTarefas);
	printf("Tempo médio por tarefa: %.6f (%.6f) ms\n", mediaTempo, desvioPadraoTempoGasto);
	printf("Fila estava vazia: %d vezes\n", contadorFilaVazia);
	return 0;
}
