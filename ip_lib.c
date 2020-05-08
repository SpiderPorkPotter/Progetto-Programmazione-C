/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

/*///////////////////FUNZIONI FEDERICO////////////////////*/

 void ip_mat_init_random(ip_mat * t, float mean, float var)
 {
	 /*k è canale(livello), h = righe(rig), w = colonne(col)*/
	 unsigned int liv, rig, col;
     float valore, risNormRand;

     for(liv = 0; liv < t -> k; liv++)
     {
         for(rig = 0; rig < t -> h; rig++)
         {
             for(col = 0; col < t -> w; col++)
             {
                 /*ottengo il valore*/
                 risNormRand = get_normal_random();
                 valore = mean + (var * risNormRand);
                 t -> data[liv][rig][col] = valore ;
                 /*per stampare il valore v
                  * printf("pro, rig, col: [%u][%u][%u] = v: %f \n", pro, rig, col, v);*/
             }
         }
     }
 }

void compute_stats(ip_mat * t)
{
	unsigned int liv, rig, col; /*k è canale(livello), h = righe(rig), w = colonne(col)*/
	int nElem;
	float totale, minimo, maximo, mid;

	/*per ogni livello. liv == 0 è stat[0]*/
	for(liv = 0; liv < t -> k; liv++)
	{
        nElem = 0;
        totale = 0;
        minimo = t -> data[0][0][0];
        maximo = t -> data[0][0][0];
		for(rig = 0; rig < t -> h; rig++)
		{
			for(col = 0; col < t -> w; col++)
			{
				/*calcolo del massimo*/
				if(t -> data[liv][rig][col] > maximo)
					maximo = t -> data[liv][rig][col];

				/*calcolo del minimo*/
				if(t -> data[liv][rig][col] < minimo)
					minimo = t -> data[liv][rig][col];

				/*calcoli necessari alla media*/
				totale += t -> data[liv][rig][col];
				nElem++;
			}
		}
		/*terminato un canale (livello), si assegnano i valori nella corrispettiva struct in stat*/
		mid = totale / nElem;
		t -> stat[liv].min = minimo;
		t -> stat[liv].max = maximo;
		t -> stat[liv].mean = mid;
	}
	/*stat da inizializzare con i valori min, max e mean (media)
	 * presenti nei vari livelli della matrice data.
	 * Quindi, dei k livelli che formano lo
	 * spessore della nostra matrice, bisogna trovare il massimo,
	 * il minimo e farne la media per ogni livello.
	 * Quindi:
	 * 		liv 1 ----> min = valore x
	 * 					max = valore xx
	 * 					media = valore xxx
	 *
	 * 		liv 2 ----> min = valore y
	 * 					max = valore yy
	 * 					media = valore zzz
	 *
	 * 		liv 3 ----> min = valore z
	 * 					max = valore zz
	 * 					media = valore zzz
	 * 		*/
}

void ip_mat_free(ip_mat *a)
{
	/*liberare data (la matrice), stat (il vettore), e poi tutta la struttura*/
	unsigned int i,j;

	/*la deallocazione parte dall'interno*/
	for (i = 0; i < a -> k; i++)
	{
		for (j = 0; j < a -> h; j++)
		{
			free(a -> data[i][j]);
		}
		free(a -> data[i]);
	}
	free(a -> data);

    /*libero ip_mat (alle 4 malloc create in ip_mat_create corrispondono le 4 free) */
    free(a);
	/*libero il vettore stat*/
	free(a -> stat);
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v)
{
	unsigned int rig, col, pro; /*indici*/

	/*variabili per capirsi meglio*/
	unsigned int righe = h;
	unsigned int colonne = w;
	unsigned int spessore = k;

	ip_mat *matrix =(ip_mat*) malloc(sizeof(ip_mat));

	/*se l'allocazione non va a buon fine*/
	if(!matrix)
		return NULL;

	/*perché ne creo una di struttura di tipo ip_mat, al cui interno va tutto*/
	/*matrix = (ip_mat) malloc(sizeof(ip_mat)*1); */
	matrix -> w = colonne;
	matrix -> h = righe;
	matrix -> k = spessore;

	float ***array;

	/*primo "cubetto"*/
	if((array = (float ***) malloc (sizeof(float ***) * spessore)))
	{
		for(pro = 0; pro < spessore; pro++)	/*ciclo per lo spessore*/
		{
			/*riga completa che parte dalla posizione rig*/
			if((array[rig] = (float **) malloc(sizeof(float*) * righe)))
			{
				/*da ogni posizione rig vengono allocate le colonne*/
				for(rig = 0; rig < righe; rig++ )
				{
					if((array[pro][rig] = (float *) malloc(sizeof(float) * colonne))){
					}
					else
						return NULL;
				}
			}
			else
				return NULL;
		}
	}
	else
		return NULL;

	/*inizializzazione col parametro v come descritto in ip_lib.h*/
	for(pro = 0; pro < spessore; pro++)
	{
		for(rig = 0; rig < righe; rig++)
		{
			for(col = 0; col < colonne; col++)
			{
				array[pro][rig][col] = v;
				/*per stampare il valore v
				 * printf("pro, rig, col: [%u][%u][%u] = v: %f \n", pro, rig, col, v);*/
			}
		}
	}
	/*perché tutta la matrice è contenuta in matrix.data*/
	matrix -> data = array;
	matrix -> stat = (stats *) malloc (sizeof(stats) * spessore);

	return matrix;
}
/*///////////////////FUNZIONI FEDERICO////////////////////*/


void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

}
