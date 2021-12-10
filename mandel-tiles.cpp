#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <queue>
#include <vector>
#include <math.h>
#include <time.h>

// Numero maximo de iteracoes maximo usado no calculo da funcao
#define MAXITER 32768
// Numero de threads trabalhadoras
#define NUMERO_THREADS_TRABALHADORAS 1

// Parametros para cada chamada a funcao de fractal
// Esses parametros indicam uma tarefa de calculo da funcao
typedef struct {
	int left; int low;  // Canto inferior esquerdo da tela (informacao lida da entrada, mas nao usada ja que nada e exibido)
	int ires; int jres; // Resolucao em pixels da area a ser computada
	double xmin; double ymin; // Canto inferior esquerdo no dominio (x, y)
	double xmax; double ymax; // Canto superior direito no dominio (x, y)
} fractal_param_t;

FILE* entrada; // Arquivo que contem a lista de blocos a serem calculados
std::queue<fractal_param_t> tarefas; // Fila de tarefas a serem consumidas pelas threads trabalhadoras
pthread_mutex_t acessoVariaveisEstatistica; // Mutex para controlar o acesso as variaveis usadas pelos calculos das estatisticas
pthread_mutex_t acessoFilaTarefas; // Mutex para controlar o acesso a fila de tarefas
pthread_cond_t preencherTarefas; // Variavel de condicao que controla se novas tarefas devem ser colocadas na fila
pthread_cond_t esperarTarefas; // Variavel de condicao que controla se novas tarefas foram colocadas na fila
int totalTarefasExecutadas = 0; // Indica o numero de tarefas executadas pelas threads trabalhadoras
int tarefasExecutadasThread[NUMERO_THREADS_TRABALHADORAS]; // Indica o numero de tarefas executadas por cada thread trabalhadora
double tempoTotalGasto = 0.0; // Indica o tempo total gasto executando as tarefas em milissegundos
std::vector<double> tempoGastoTarefas; // Indica o tempo gasto executando cada tarefa em milissegundos
int contadorFilaVazia = 0; // Indica o numero de vezes que uma thread trabalhadora encontrou a fila de tarefas vazia

// Realiza os cálculos indicados pela tarefa recebida por parâmetro
void realizarCalculo(fractal_param_t* p) {
	// Declara as variaveis necessarias para os calculos da tarefa
	double dx, dy;
	double x, y, u, v, u2, v2;
	
	// Calcula o valor de dx e dy
	dx = ( p->xmax - p->xmin ) / p->ires;
	dy = ( p->ymax - p->ymin ) / p->jres;
	
	// Escaneia cada ponto na area retangular da tarefa
	// Cada ponto representa um numero complexo (x + yi)
	// Itera por esse numero complexo
	for(int j = 0; j < p->jres; j++) {
		for(int i = 0; i <= p->ires; i++) {
			x = i * dx + p->xmin; // c_real
			u = u2 = 0; // z_real
			y = j * dy + p->ymin; // c_imaginario
			v = v2 = 0; // z_imaginario

			// Calcula se c(c_real + c_imaginario) pertence ao conjunto de Mandelbrot e desenha um pixel na coordenada (i, j) de acordo
			// Se foi alcançado o numero maximo de iteracoes e se a distancia da origem e maior do que dois, encerra o loop
			for(int k=0; (k < MAXITER) && ((u2+v2) < 4); ++k) {
				// Calcula a funcao de Mandelbrot
				// z = z*z + c em que z e um numero complexo

				// imag = 2*z_real*z_imaginary + c_imaginary
				v = 2 * u * v + y;
				// real = z_real^2 - z_imaginary^2 + c_real
				u  = u2 - v2 + x;
				u2 = u * u;
				v2 = v * v;
			}
		}
	}
}

// Thread responsavel por consumir tarefas da fila de tarefas e executa-las, gerando dados para as estatisticas de execucao
void* threadTrabalhadora(void* param) {
	// Extrai o id da thread do parametro recebido pela thread
	long threadId = (long) param;
	// Indica se a thread leu uma flag de end of work
	bool leuEOW = false;
	// Executa tarefas enquanto nao for encontrada uma tarefa que indica end of work
	while(!leuEOW) {
		// Declara flag que indica se a thread encontrou a fila vazia
		bool filaFicouVazia = false;
		// Obtem o mutex de acesso a fila de tarefas
		pthread_mutex_lock(&acessoFilaTarefas);
		// Enquanto a lista de tarefas estiver vazia
		while(tarefas.empty()) {
			// Indica que a thread encontrou a lista de tarefas vazia
			filaFicouVazia = true;
			// Espera a thread leitora preencher a fila de tarefas
			pthread_cond_wait(&esperarTarefas, &acessoFilaTarefas);
		}
		// Obtem a proxima tarefa da fila
		fractal_param_t p =	tarefas.front();
		// Remove essa tarefa da fila
		tarefas.pop();
		// Caso a fila de tarefas possua o mesmo numero de tarefas que de threads trabalhadoras
		if(tarefas.size() == NUMERO_THREADS_TRABALHADORAS) {
			// Indica para a thread leitora que ela deve preencher a fila de tarefas
			pthread_cond_signal(&preencherTarefas);
		}
		// Libera o mutex de acesso a fila de tarefas
		pthread_mutex_unlock(&acessoFilaTarefas);
		// Caso a tarefa seja indicativo de end of work
		if(p.jres == 0 && p.ires == 0) {
			// Indica que a thread leu um end of work
			leuEOW = true;
		}
		// Caso a tarefa nao seja de end of work
		if(!leuEOW) {
			// Inicia a contagem do tempo gasto para executar da tarefa
			clock_t inicio = clock();
			// Realiza os cálculos da tarefa
			realizarCalculo(&p);
			// Finaliza a contagem do tempo gasto para executar a tarefa
			clock_t fim = clock();
			// Calcula o valor do tempo gasto em milissegundos
			double tempoGasto = ((double)fim - inicio) / CLOCKS_PER_SEC * 1000;
			// Obtem o mutex que controla o acesso as variaveis estatisticas
			pthread_mutex_lock(&acessoVariaveisEstatistica);
			// Salva que mais uma tarefa foi executada por essa thread trabalhadora
			totalTarefasExecutadas+=1;
			tarefasExecutadasThread[threadId]++;
			// Salva o tempo gasto para executar essa tarefa
			tempoTotalGasto += tempoGasto;
			tempoGastoTarefas.push_back(tempoGasto);
			// Caso essa thread tenha encontrado a fila vazia
			if(filaFicouVazia) {
				// Salva que essa thread encontrou a fila vazia ao executar essa tarefa
				contadorFilaVazia++;
			}
			// Libera o mutex que controla o acesso as variaveis estatisticas
			pthread_mutex_unlock(&acessoVariaveisEstatistica);
		}
	}
	// Encerra a execucao da thread leitora
	pthread_exit(NULL);
}

// Thread responsavel por ler as tarefas do arquivo e por preencher a fila de tarefas
void* threadLeitora(void* param) {
	// Declara uma variavel do tipo que representa os parametros de uma tarefa
	fractal_param_t p;
	// Obtem o mutex do acesso a fila de tarefas
	pthread_mutex_lock(&acessoFilaTarefas);
	// Flag para indicar se todo o arquivo ja foi lido
	bool acabouArquivo = false;
	// Le os quatro primeiros parametros de uma tarefa do arquivo
	int n = fscanf(entrada,"%d %d %d %d",&(p.left),&(p.low),&(p.ires),&(p.jres));
	// Se o arquivo nao possui mais dados, indica que ele foi totalmente lido
	if (n == EOF) acabouArquivo = true;
	while(!acabouArquivo) {
		// Verifica se ocorreu um erro no fscanf
		if(n != 4) {
			// Indica que se teve um erro no fscanf e encerra o programa
			perror("fscanf(left,low,ires,jres)");
			exit(-1);
		}
		// Le os quatro ultimos parametros de uma tarefa do arquivo
		n = fscanf(entrada,"%lf %lf %lf %lf", &(p.xmin),&(p.ymin),&(p.xmax),&(p.ymax));
		// Verifica se ocorreu um erro no fscanf
		if(n != 4) {
			// Indica que se teve um erro no fscanf e encerra o programa
			perror("fcanf(xmin,ymin,xmax,ymax)");
			exit(-1);
		}
		// Adiciona a tarefa lida na fila de tarefas
		tarefas.push(p);
		// Verifica se a fila já possui 4 vezes a quantidade de threads de tarefas
		if(tarefas.size() >= 4*NUMERO_THREADS_TRABALHADORAS) {
			// Caso a fila ja esteja cheia, indica que a fila foi preenchida caso alguma thread esteja esperando por tarefas
			pthread_cond_broadcast(&esperarTarefas);
			// Apos isso, espera algum thread trabalhadora indicar que a fila de tarefas deve ser preenchida novamente
			pthread_cond_wait(&preencherTarefas, &acessoFilaTarefas);
		}
		// Le os quatro primeiros parametros de uma tarefa do arquivo
		n = fscanf(entrada,"%d %d %d %d",&(p.left),&(p.low),&(p.ires),&(p.jres));
		// Se o arquivo nao possui mais dados, indica que ele foi totalmente lido
		if(n == EOF) acabouArquivo = true;
	}
	// Salva na variavel dos parametros a indicacao de que essa tarefa indica um end of work para as threads trabalhadoras
	p.ires = 0;
	p.jres = 0;
	// Coloca na fila uma tarefa de end of work para cada thread trabalhadora
	for(int i=0; i < NUMERO_THREADS_TRABALHADORAS; i++) {
		tarefas.push(p);
	}
	// Libera o mutex de acesso a fila de tarefas
	pthread_mutex_unlock(&acessoFilaTarefas);
	// Encerra a execucao da thread leitora
	pthread_exit(NULL);
}

int main(int argc, char* argv[]) {

	// Verifica se o programa recebeu o arquivo de entrada como argumento
	if(argc!=2) {
		// Imprime o uso correto do programa e encerra sua execucao
		fprintf(stderr,"usage %s filename\n", argv[0] );
		exit(-1);
	}
	
	// Verifica se o arquivo recebido por argumento pode ser aberto
	if((entrada=fopen(argv[1],"r"))==NULL) {
		// Caso ocorra um erro ao abrir o arquivo, encerra o programa
		perror("fopen");
		exit(-1);
	}

	// Inicializa os dois mutex usados para sincronizacao
	pthread_mutex_init(&acessoVariaveisEstatistica, NULL);
	pthread_mutex_init(&acessoFilaTarefas, NULL);
	// Inicializa as duas variaveis de condicao usadas para sincronizacao
	pthread_cond_init(&preencherTarefas, NULL);
	pthread_cond_init(&esperarTarefas, NULL);

	// Cria a thread responsavel por ler o arquivo e popular a fila de tarefas
	pthread_t threadLeitura;
	pthread_create(&threadLeitura, NULL, threadLeitora, NULL);

	// Cria as threads responsaveis por consumir as tarefas da fila e executa-las
	pthread_t threadsTrabalhadoras[NUMERO_THREADS_TRABALHADORAS];
	for(long i=0; i<NUMERO_THREADS_TRABALHADORAS; i++) {
		pthread_create(&threadsTrabalhadoras[i], NULL, threadTrabalhadora, (void*) i);
	}

	clock_t inicio = clock();

	// Espera a finalizacao da thread de leitura
	pthread_join(threadLeitura, NULL);

	// Espera a finalizacao de todas as threads trabalhadoras
	for(int i=0; i<NUMERO_THREADS_TRABALHADORAS; i++) {
		pthread_join(threadsTrabalhadoras[i], NULL);
	}

	clock_t fim = clock();

	double tempoGasto = ((double)fim - inicio) / CLOCKS_PER_SEC;
	printf("Tempo total gasto = %.6f\n", tempoGasto);

	// Calcula as estasticas da execucao das tarefas
	// Calcula a media do numero de tarefas executadas por cada thread
	double mediaTarefas = (double) totalTarefasExecutadas/(double) NUMERO_THREADS_TRABALHADORAS;

	// Calcula o desvio padrao do numero de tarefas executadas por cada thread
	double desvioPadraoTarefas = 0;
	for(int i=0; i<NUMERO_THREADS_TRABALHADORAS; i++) {
		desvioPadraoTarefas += (tarefasExecutadasThread[i] - mediaTarefas)*(tarefasExecutadasThread[i] - mediaTarefas);
	}
	desvioPadraoTarefas /= (double) NUMERO_THREADS_TRABALHADORAS;
	desvioPadraoTarefas = sqrt(desvioPadraoTarefas);

	// Calcula a media de tempo gasto para executar cada tarefa
	double mediaTempo = tempoTotalGasto / (double) totalTarefasExecutadas;

	// Desvio padrão do tempo gasto para executar cada tarefa
	double desvioPadraoTempoGasto = 0;
	for(int i=0; i<totalTarefasExecutadas; i++) {
		desvioPadraoTempoGasto += (tempoGastoTarefas[i] - mediaTempo)*(tempoGastoTarefas[i] - mediaTempo);
	}
	desvioPadraoTempoGasto /= (double) totalTarefasExecutadas;
	desvioPadraoTempoGasto = sqrt(desvioPadraoTempoGasto);

	// Imprime as estatisticas geradas pelas threads
	printf("Tarefas: total = %d;  média por trabalhador = %lf(%lf)\n", totalTarefasExecutadas, mediaTarefas, desvioPadraoTarefas);
	printf("Tempo médio por tarefa: %.6f (%.6f) ms\n", mediaTempo, desvioPadraoTempoGasto);
	printf("Fila estava vazia: %d vezes\n", contadorFilaVazia);

	// Destroi os dois mutex usados para sincronizacao
	pthread_mutex_destroy(&acessoVariaveisEstatistica);
	pthread_mutex_destroy(&acessoFilaTarefas);
	// Destroi as duas variaveis de condicao usadas para sincronizacao
	pthread_cond_destroy(&preencherTarefas);
	pthread_cond_destroy(&esperarTarefas);

	return 0;
}
