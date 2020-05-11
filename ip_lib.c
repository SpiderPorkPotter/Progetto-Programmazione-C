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

/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/
/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
ip_mat * ip_mat_to_gray_scale(ip_mat * in)
{
    /*k è canale(livello), h = righe(rig), w = colonne(col)*/
    unsigned int liv, rig, col;
    float totale, media;
    /*ne creo una nuova con tutti gli elementi a 0*/
    ip_mat *out = ip_mat_create(in -> h, in -> w, in -> k, 0);

    for(rig = 0; rig < in -> w; rig++)
    {
        for(col = 0; col < in -> w; col++)
        {
            /*calcolo la media di ogni pixel e la scrivo sui pixel del primo canale*/
            for(liv = 0; liv < in -> k; liv++)
            {
                totale += in -> data[liv][rig][col];
            }
            /*finito il ciclo, è stata calcolata la somma dello stesso pixel sui livelli
            quindi ora non resta che calcolare la media e assegnarla alla cella in cui stiamo "puntando"
            in questo momento. Così poi passiamo alla successiva e ripetiamo*/
            media = totale / in -> k;
            out -> data[0][rig][col] = media;
        }
    }
    /*finito questo, significa che abbiamo iterato per tutto il primo canale, assegnando ad ogni
    pixel la sua media. Ora bisogna assegnare i corrispettivi valori dei pixel anche ai livelli successivi*/

    for(rig = 0; rig < in -> w; rig++)
    {
        for(col = 0; col < in -> w; col++)
        {
            /*ai canali successivi al primo (lo 0) assegno i valori presenti nel primo. Ovvero, copio i valori
            del primo livello nel secondo e nel terzo, cella per cella, scorrendo i canali. Poi passo all'elemento
            successivo del primo livello (col++) e ripeto*/
            for(liv = 1; liv < in -> k; liv++)
            {
                out -> data[liv][rig][col]= in -> data[0][rig][col];
            }
        }
    }
    /*per riempire il vettore stats*/
    compute_stats(out);
    /*l'immagine fornita in ingresso nella funzione non serve più, quindi può essere liberata*/
    ip_mat_free(in);
    return out;
}

/*///////////////////FUNZIONI FEDERICO////////////////////*/


/* ###########################
   Inizio funzioni parte 1 TOM
   ########################### */

/* IP MAT COPY - crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in) {
    ip_mat *x;  /* ho creato il puntatore ad una nuova struttura ip_mat */
    unsigned int liv, rig, col;
    /* ho richiamato la funzione create per assegnare le tre dimensioni ed allocare stats e data */
    x = ip_mat_create(in->h, in->w, in->k, 1.0);

    /* ora assegnamo i valori di data di in a x */
    for ( liv=0; liv<in->k; liv++ ) {
        for ( rig=0; rig<in->h; rig++ ) {
            for ( col=0; col<in->w; col++ ) {
                x->data[liv][rig][col] = in->data[liv][rig][col];
            }
        }
    }

    /* ora assegnamo i valori di stats di in a x */
    for ( liv=0; liv<in->k; liv++ ) {
        x->stat[liv].min = in->stat[liv].min;
        x->stat[liv].max = in->stat[liv].max;
        x->stat[liv].mean = in->stat[liv].mean;
    }

    return x;
}


 ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
     ip_mat *x;  /* ho creato il puntatore ad una nuova struttura ip_mat */
     unsigned int liv, rig, col;
     x = NULL;

     /* controlla se:
        col e row end sono maggiori del loro corrispettivo starts
        se entrambi gli start sono maggiori o uguali a zero
        se entrambi gli end sono massimo grandi quanto quelli presenti in t
     */
     if ( (row_start <= row_end)  &&  (col_start <= col_end)  &&  (row_start >= 0)  &&  (col_start >= 0)  &&  (row_end <= t->h)  &&  (col_end <= t->w) ){

         x = ip_mat_create( (row_end - row_start) , (col_end - col_start) ,t->k, 1.0); /* istanzio la nuova truttura ip_mat con il numero di righe e colonne richieste */

         /* ora vado a riempire le righe e colonne e k di data*/
         for ( liv=0; liv<t->k; liv++ ) {
             for ( rig=row_start; rig<row_end; rig++ ) {
                 for ( col=col_start; col<col_end; col++ ) {
                     x->data[liv][rig-row_start][col-col_start] = t->data[liv][rig][col];
                 }
             }
         }

         /* faccio elaborare gli stats */
         compute_stats(x);
     }


     /* se non si entra nell'if vuol dire che c'è un problema e quindi restituisce NULL, altrimenti ritorna il puntatore del nuovo ip_mat */
     return x;
 }


  ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione) {
      ip_mat *x; /* ho creato il puntatore ad una nuova struttura ip_mat */
      unsigned int liv, rig, col;
      x = NULL;


      /* concateno l'altezza */
      if ( dimensione == 0 ) {
          if ( (a->w == b->w)  &&  (a->k == b->k) ){
            ip_mat_create( (a->h + b->h), a->w, a->k, 0.0 );
            /* ora vado a riempire le righe e colonne e k di data*/
                for ( liv=0; liv<x->k; liv++ ) {
                    for ( rig=0; rig<x->h; rig++ ) {
                        for ( col=0; col<x->w; col++ ) {
                            if ( rig < a->h ) {
                                x->data[liv][rig][col] = a->data[liv][rig][col];
                            }
                            else {
                                x->data[liv][rig][col] = b->data[liv][rig-(a->h)][col];
                            }
                        }
                    }
                }
            }
      }

      /* Concateno la larghezza */
      if ( dimensione == 1 ) {
          if ( (a->h == b->h)  &&  (a->k == b->k) ){
            ip_mat_create( a->h, (a->w + b->w), a->k, 0.0 );
            /* ora vado a riempire le righe e colonne e k di data*/
                for ( liv=0; liv<x->k; liv++ ) {
                    for ( rig=0; rig<x->h; rig++ ) {
                        for ( col=0; col<x->w; col++ ) {
                            if ( col < a->w ) {
                                x->data[liv][rig][col] = a->data[liv][rig][col];
                            }
                            else {
                                x->data[liv][rig][col] = b->data[liv][rig][col - (a->w)];
                            }
                        }
                    }
                }
            }
      }

      /* Concateno la profondità (k) */
      if ( dimensione == 2 ) {
          if ( (a->h == b->h)  &&  (a->k == b->k) ){
            ip_mat_create( a->h, a->w, (a->k + b->k), 0.0 );
            /* ora vado a riempire le righe e colonne e k di data*/
                for ( liv=0; liv<x->k; liv++ ) {
                    for ( rig=0; rig<x->h; rig++ ) {
                        for ( col=0; col<x->w; col++ ) {
                            if ( liv < a->k ) {
                                x->data[liv][rig][col] = a->data[liv][rig][col];
                            }
                            else {
                                x->data[liv][rig][col] = b->data[liv - a->k][rig][col];
                            }
                        }
                    }
                }
            }
      }

      return x;
  }

/* ###########################
   Fine funzioni parte 1 TOM
   ########################### */


/*//////////////Yamina & Christian FUNCTIONS//////////////*/
/***     PARTE 1: OPERAZIONE MATEMATICHE FRA IP_MAT     ***/

/*Sum of two ip_mat (all dimensions must be the same)*/
ip_mat *ip_mat_sum (ip_mat *a, ip_mat *b)
{
    unsigned int liv, rig, col;
    ip_mat *out; /*create the ip_mat output*/

    out = NULL; /*initialize the new ip_mat*/

    if(a->h == b->h && a->w == b->w && a->k == b->k) /*control if all dimensions are equal*/
    {
        for(liv = 0; liv < a->k; liv++) /*cycle to scroll the depth*/
            for(rig = 0; rig < a->h; rig++) /*cycle to scroll the height*/
                for(col = 0; col < a->w; col++) /*cycle to scroll the width*/
                    out->data[liv][rig][col] = a->data[liv][rig][col] + b->data[liv][rig][col];
    }
    return out;
}

/*Subtraction of two ip_mat (all dimensions must be the same)*/
ip_mat *ip_mat_sub (ip_mat *a, ip_mat *b)
{
    unsigned int liv, rig, col;
    ip_mat *out; /*create the ip_mat output*/

    out = NULL; /*initialize the new ip_mat*/

    if(a->h == b->h && a->w == b->w && a->k == b->k) /*control if all dimensions are equal*/
    {
        for(liv = 0; liv < a->k; liv++) /*cycle to scroll the depth*/
            for(rig = 0; rig < a->h; rig++) /*cycle to scroll the height*/
                for(col = 0; col < a->w; col++) /*cycle to scroll the width*/
                    out->data[liv][rig][col] = a->data[liv][rig][col] - b->data[liv][rig][col];
    }
    return out;
}

/*Multiply a scalar to a ip_mat*/
ip_mat * ip_mat_mul_scalar (ip_mat *a, float c)
{
    unsigned int liv, rig, col;
    ip_mat *out; /*create the ip_mat output*/

    out = NULL; /*initialize the new ip_mat*/

    for(liv = 0; liv < a->k; liv++) /*cycle to scroll the depth*/
        for(rig = 0; rig < a->h; rig++) /*cycle to scroll the height*/
            for(col = 0; col < a->w; col++) /*cycle to scroll the width*/
                out->data[liv][rig][col] = a->data[liv][rig][col] * c;

    return out;
}

/*Add a scalar to a ip_mat*/
ip_mat * ip_mat_add_scalar (ip_mat *a, float c)
{
    unsigned int liv, rig, col;
    ip_mat *out; /*create the ip_mat output*/

    out = NULL; /*initialize the new ip_mat*/

    for(liv = 0; liv < a->k; liv++) /*cycle to scroll the depth*/
        for(rig = 0; rig < a->h; rig++) /*cycle to scroll the height*/
            for(col = 0; col < a->w; col++) /*cycle to scroll the width*/
                out->data[liv][rig][col] = a->data[liv][rig][col] + c;

    return out;
}

/*Mean between two ip_mat*/
/*    ????? Secondo voi è meglio usare la funzione ip_mat_sum ????	*/
ip_mat * ip_mat_mean (ip_mat *a, ip_mat *b)
{
    unsigned int liv, rig, col;
    ip_mat *out; /*create the ip_mat output*/

    out = NULL; /*initialize the new ip_mat*/

    if(a->h == b->h && a->w == b->w && a->k == b->k) /*control if all dimensions are equal*/
    {
        for(liv = 0; liv < a->k; liv++) /*cycle to scroll the depth*/
            for(rig = 0; rig < a->h; rig++) /*cycle to scroll the height*/
                for(col = 0; col < a->w; col++) /*cycle to scroll the width*/
                    out->data[liv][rig][col] = (a->data[liv][rig][col] + b->data[liv][rig][col])/2;
    }
    return out;
}

/*//////////////Yamina & Christian FUNCTIONS//////////////*/

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
