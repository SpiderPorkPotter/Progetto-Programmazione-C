/*
 Created by Sebastiano Vascon on 23/03/20.

 WLF
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

/*///////////////////FUNZIONI FEDERICO PARTE 1////////////////////*/

 void ip_mat_init_random(ip_mat * t, float mean, float var)
 {
	 /*k è canale(livello), h = righe(rig), w = colonne(col)*/
	 unsigned int liv, rig, col;
     float risNormRand;

     if(t == NULL)
     {
         printf("Errore! IPMAT non valida/e\n");
         exit(1);
     }

     for(liv = 0; liv < t -> k; liv++)
     {
         for(rig = 0; rig < t -> h; rig++)
         {
             for(col = 0; col < t -> w; col++)
             {
                 /*ottengo il valore*/
                 risNormRand = get_normal_random(mean, var);
                 t -> data[liv][rig][col] = risNormRand ;
             }
         }
     }
     compute_stats(t);
 }

void compute_stats(ip_mat * t)
{
	unsigned int liv, rig, col; /*k è canale(livello), h = righe(rig), w = colonne(col)*/
	int nElem;
	float totale, minimo, maximo, mid;

    if(t == NULL)
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

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
	/*liberare data (la matrice), stat (il vettore), e tutta la struttura*/
	unsigned int liv,rig;
    if(a != NULL) {
        /*libera stat*/
        if(a -> stat != NULL)
        free(a -> stat);
        /*elimina prima per ogni cella delle righe, la colonna.
        Terminato il ciclo interno elimina tutte le righe e quindi la matrice.
        Passa al livello (canale) successivo e ripete*/
        if(a -> data != NULL)
        {
            for (liv = 0; liv < a -> k; liv++)
            {
                for (rig = 0; rig < a -> h; rig++)
                {
                    /*if(a -> data[liv][rig] != NULL)*/
                        free(a -> data[liv][rig]);
                }
                if(a -> data[liv] != NULL)
                    free(a -> data[liv]);
            }
            free(a -> data);
        }
        free(a);
    }
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v)
{
    unsigned int liv, rig, col; /*indici*/
    ip_mat* Ipmat;

    /*se l'allocazione non va a buon fine*/
    Ipmat = (ip_mat*)malloc(sizeof(ip_mat));

    if(Ipmat == NULL)
    {
        printf("Errore nell'allocazione struttura ipmat\n" );
        exit(1);
    }

    /*perché ne creo una di struttura di tipo ip_mat, al cui interno va tutto*/
    Ipmat -> w = w; /*colonne*/
    Ipmat -> h = h; /*righe*/
    Ipmat -> k = k; /*livello (o canale)*/

    /*se l'allocazione dei tre livelli (canali) va a buon fine*/
    if( ( Ipmat -> data = (float***) malloc(k * sizeof(float**)) ) )
    {

        for(liv = 0; liv < k; liv++ ) /*per ogni livello*/
        {
            if( ( Ipmat -> data[liv] = (float**) malloc(h * sizeof(float*)) ) ) /*alloca h righe (puntatori)*/
            {
                for(rig = 0; rig < h; rig++)    /*per ogni puntatore*/
                {
                    if( ! ( Ipmat -> data[liv][rig] = (float*) malloc(w * sizeof(float)) ) )/*alloca una colonna composta da float*/
                    {
                        printf("Errore! Impossibile allocare!\n");
                        free(Ipmat);
                        exit(1);
                    }
                }
            }
            else
            {
                printf("Errore! Impossibile allocare!\n");
                free(Ipmat);
                exit(1);
            }
        }
    }
    else
    {
        printf("Errore! Impossibile allocare!\n");
        free(Ipmat);
        exit(1);
    }

    /*inizializzazione col parametro v come descritto in ip_lib.h*/
    for(liv = 0; liv < k; liv++)
    {
        for(rig = 0; rig < h; rig++)
        {
            for(col = 0; col < w; col++)
            {
                Ipmat -> data[liv][rig][col] = v;
            }
        }
    }
    /*perché tutta la matrice è contenuta in matrix.data*/
    /*Ipmat -> data = matrix;*/
    Ipmat -> stat = (stats*) malloc(k * sizeof(stats));
    if(Ipmat -> stat == NULL)
    {
        printf("Errore nell'allocazione struttura stat\n");
        exit(1);
    }
    /* inizializzazione valori vettore stats */
    for ( liv=0; liv < k; liv++ )
    {
        Ipmat->stat[liv].min = v;
        Ipmat->stat[liv].max = v;
        Ipmat->stat[liv].mean = v;
    }

    return Ipmat;
}

/*///////////////////FUNZIONI FEDERICO PARTE 1////////////////////*/


/* ###########################
   Inizio funzioni parte 1 TOM
   ########################### */

/* IP MAT COPY - crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in) {
    ip_mat *x;  /* ho creato il puntatore ad una nuova struttura ip_mat */
    unsigned int liv, rig, col;

    if(in == NULL)
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

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
     unsigned int liv, rig, col, zero=0;
     x = NULL;

     if(t == NULL)
     {
         printf("Errore! IPMAT non valida/e\n");
         exit(1);
     }

     /* controlla se:
        col e row end sono maggiori del loro corrispettivo starts
        se entrambi gli start sono maggiori o uguali a zero
        se entrambi gli end sono massimo grandi quanto quelli presenti in t
     */
     if ( (row_start <= row_end)  &&  (col_start <= col_end)  &&  (row_start >= zero)  &&  (col_start >= zero)  &&  (row_end <= t->h)  &&  (col_end <= t->w) ){

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

      if((a == NULL || b == NULL) || (a == NULL && b == NULL) ||  ( dimensione < 0)  ||  ( dimensione > 2 ) )
      {
          printf("Errore! IPMAT non valida/e\n");
          exit(1);
      }
      x = NULL;

      /* concateno l'altezza */
      if ( dimensione == 0 ) {
          if ( (a->w == b->w)  &&  (a->k == b->k) ){
            x = ip_mat_create( (a->h + b->h), a->w, a->k, 0.0 );
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
            x = ip_mat_create( a->h, (a->w + b->w), a->k, 0.0 );
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
            x = ip_mat_create( a->h, a->w, (a->k + b->k), 0.0 );
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
    ip_mat *out;

    if((a == NULL || b == NULL) || (a == NULL && b == NULL))
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

    out = ip_mat_create(a->h, a->w, a->k, 0);  /*create the ip_mat output*/

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
    ip_mat *out;

    if((a == NULL || b == NULL) || (a == NULL && b == NULL))
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

    out = ip_mat_create(a->h, a->w, a->k, 0);  /*create the ip_mat output*/

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
    ip_mat *out;

    if(a == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

    out = ip_mat_create(a->h, a->w, a->k, 0);  /*create the ip_mat output*/

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
    ip_mat *out;

    if(a == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

    out = ip_mat_create(a->h, a->w, a->k, 0);  /*create the ip_mat output*/

    for(liv = 0; liv < out->k; liv++) /*cycle to scroll the depth*/
    {
        for(rig = 0; rig < out->h; rig++) /*cycle to scroll the height*/
        {
            for(col = 0; col < out->w; col++) /*cycle to scroll the width*/
            {
                out->data[liv][rig][col] = (a->data[liv][rig][col]) + c;
            }
        }
    }
    clamp(out,0,255);
    return out;
}

/*Mean between two ip_mat*/
/*    ????? Secondo voi è meglio usare la funzione ip_mat_sum ????	*/
ip_mat * ip_mat_mean (ip_mat *a, ip_mat *b)
{
    unsigned int liv, rig, col;
    ip_mat *out;

    if((a == NULL || b == NULL) || (a == NULL && b == NULL))
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

    out = ip_mat_create(a->h, a->w, a->k, 0); /*create the ip_mat output*/

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

/* ###############
   ### PARTE 2 ###
   ############### */

/*
IP_MAT_TO_GRAY_SCALE
Autore: Federico
Descrizione: prende un'immagine colorata e la rende in scala di grigi. Si calcola la media
            dei tre canali per ogni pixel e viene assegnata al primo livello per ogni pixel.
            Calcolato per tutto il primo livello, si copiano i valori anche sugli altri 2 livelli.
*/

ip_mat * ip_mat_to_gray_scale(ip_mat * in)
{
    /*k è canale(livello), h = righe(rig), w = colonne(col)*/
    unsigned int liv, rig, col;
    float totale = 0, media = 0;
    ip_mat *out;

    if(in == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

    /*ne creo una nuova con tutti gli elementi a 0*/
    out = ip_mat_create(in -> h, in -> w, in -> k, 0);

    for(rig = 0; rig < in -> h; rig++)
    {
        for(col = 0; col < in -> w; col++)
        {
            /*calcolo la media di ogni pixel e la scrivo sui pixel del primo canale*/
            totale = 0;
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

    for(rig = 0; rig < in -> h; rig++)
    {
        for(col = 0; col < in -> w; col++)
        {
            /*ai canali successivi al primo (lo 0) assegno i valori presenti nel primo. Ovvero, copio i valori
            del primo livello nel secondo e nel terzo, cella per cella, scorrendo i canali. Poi passo all'elemento
            successivo del primo livello (col++) e ripeto*/
            for(liv = 1; liv < in -> k; liv++)
            {
                out -> data[liv][rig][col]= out -> data[0][rig][col];
            }
        }
    }
    /*per riempire il vettore stats*/
    compute_stats(out);
    /*l'immagine fornita in ingresso nella funzione non serve più, quindi può essere liberata*/
    /*ip_mat_free(in); c'è già sul main*/
    return out;
}

/*
IP_MAT_BLEND
Autore: Tom
Descrizione: Unisce due immagini assieme, queste devono essere della stessa dimensione (anche i canali).
	     La funzione non modifica i parametri in ingresso ma crea una nuova ip_mat su cui scrive i risultati delle computazioni.
*/
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
    ip_mat *x;
    unsigned int liv, rig, col;
    x = NULL;

    /*
    if ( b == NULL )
        printf("\n\n\na == null\n\n\n");
    */
    if( ((a == NULL) || (b == NULL)) || ((a->h != b->h) || (a->w != b->w) || (a->k != b->k)) )
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

    /* ho richiamato la funzione create per assegnare le tre dimensioni ed allocare stats e data */
    x = ip_mat_create(a->h, a->w, a->k, 100);  /*prima era scritto in, ma non esiste, allora ho messo a per coerenza con l'if. Controllare se è corretto*/

    /* ora assegnamo i valori di data di in a x */
    for ( liv=0; liv<x->k; liv++ ) {
        for ( rig=0; rig<x->h; rig++ ) {
            for ( col=0; col<x->w; col++ ) {
                x->data[liv][rig][col] = (alpha * (a->data[liv][rig][col])) + ((1-alpha) * (b->data[liv][rig][col]));
            }
        }
    }

    /* calcolo le statistiche su x */
    compute_stats(x);

    return x;
}

/* IP_MAT_BRIGHTEN
Autore: Yamina
Descrizione: aumenta la luminosità dell'immagine, aggiunge ad ogni pixel un certo valore.
	     I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
	     all'interno di una nuova ip_mat.
*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright)
{
    ip_mat *out;

    if(a == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

    out = ip_mat_add_scalar(a, bright);
    compute_stats(out);
    return out;
}


/*
	IP_MAT_CORRUPT
	Autori: Fatta da tutti assieme in videochiamata
*/
ip_mat * ip_mat_corrupt(ip_mat * in, float amount) {
    ip_mat *out;  /* ho creato il puntatore ad una nuova struttura ip_mat */
    unsigned int liv, rig, col;
    float ris;

    if(in == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

    /* ho richiamato la funzione copy per copiare la ip_mat in entrata */
    out = ip_mat_copy(in);

    /* ora assegnamo i valori di data di in a x */
    for ( liv=0; liv<in->k; liv++ ) {
        for ( rig=0; rig<in->h; rig++ ) {
            for ( col=0; col<in->w; col++ ) {
                ris = (in->data[liv][rig][col]) +  get_normal_random(0, amount/2);
                out->data[liv][rig][col] = ris;
            }
        }
    }
    return out;
}

/* ####################
   ### FINE PARTE 2 ###
   #################### */

/* ###############
   ### PARTE 3 ###
   ############### */

/*
IP_MAT_PADDING
Autore: Tom
Descrizione: Funzione per aggiungere un padding ad una ip_mat passata.
	     Non lavora direttamente sulla ip_mat passata in input ma ne genera una copia (x) su cui viene elaborato il tutto.
*/
ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w) {
    ip_mat *x; /* ip_mat dove copierò a e dove eseguirò tutte le operazioni */
    unsigned int liv, rig, col;

    if(a == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }


    x = ip_mat_copy(a);

    for ( liv=0; liv<x->k; liv++ ) {
        /* sezione superiore */
        for ( rig=0; rig<pad_h; rig++ ) {
            for ( col=0; col<x->w; col++ ) {
                x->data[liv][rig][col] = 0.0;
            }
        }
        /* sezioni laterali - parte sx */
        for ( rig=pad_h; rig<((x->h)-pad_h); rig++ ) {
            for ( col=0; col<((x->w)-pad_w); col++ ) {
                x->data[liv][rig][col] = 0.0;
            }
        }
        /* sezioni laterali - parte dx */
        for ( rig=pad_h; rig<((x->h)-pad_h); rig++ ) {
            for ( col=((x->w)-pad_w); col<x->w; col++ ) {
                x->data[liv][rig][col] = 0.0;
            }
        }
        /* sezione inferiori */
        for ( rig=(x->h)-pad_h; rig<x->h; rig++ ) {
            for ( col=0; col<x->w; col++ ) {
                x->data[liv][rig][col] = 0.0;
            }
        }
    }

    return x;
}

/*
CREATE_SHARPEN_FILTER
Autore: Federico
Descrizione: Crea un filtro di sharpening
*/
ip_mat * create_sharpen_filter()
{
    unsigned int liv = 0, rig, col, dim = 3;
    ip_mat *filter;
    filter = ip_mat_create(dim, dim, dim, 0);

    if(filter == NULL)
    {
        printf("Errore! Impossibile allocare Ipmat\n");
        exit(1);
    }

    /*per lo 0:

     0   -1    0
    -1    5   -1
     0   -1    0

     metto uno 0 ogni volta che sto leggendo una cella che fa parte dei 4 angoli*/
    for(rig = 0; rig < dim; rig++)
    {
        for(col = 0; col < dim; col++)
        {
            if ( (rig == 0 && col == 0) || (rig == 0 && col == dim -1) ||
                 (rig == dim -1 && col == 0) ||  (rig == dim -1 && col == dim -1)
                )
            {
                filter -> data[liv][rig][col] = 0;
            }
        }
    }

    /*per il -1*/

    for(rig = 0; rig < dim; rig++)
    {
        for(col = 0; col < dim; col++)
        {
            if ( (rig == 0 && col == 1) || (rig == 1 && col == 0) ||
                 (rig == dim -1 && col == 1) ||  (rig == 1 && col == dim -1)
                )
            {
                filter -> data[liv][rig][col] = -1;
            }
        }
    }

    /*per il 5*/

    filter -> data[0][1][1] = 5;

    /*copia i valori anche negli altri due livelli*/
    for(rig = 0; rig < filter -> h; rig++)
    {
        for(col = 0; col < filter -> w; col++)
        {
            /*ai canali successivi al primo (lo 0) assegno i valori presenti nel primo. Ovvero, copio i valori
            del primo livello nel secondo e nel terzo, cella per cella, scorrendo i canali. Poi passo all'elemento
            successivo del primo livello (col++) e ripeto*/
            for(liv = 1; liv < filter -> k; liv++)
            {
                filter -> data[liv][rig][col] = filter -> data[0][rig][col];
            }
        }
    }
    compute_stats(filter);

    return filter;
}

/*
CREATE_EDGE_FILTER
Autore: Federico
Descrizione: Crea un filtro per rilevare i bordi
*/

ip_mat * create_edge_filter()
{
    unsigned int liv = 0, rig, col, dim = 3;
    ip_mat *filter;
    filter = ip_mat_create(dim, dim, dim, 0);

    if(filter == NULL)
    {
        printf("Errore! Impossibile allocare Ipmat\n");
        exit(1);
    }

    for(rig = 0; rig < dim; rig++)
    {
        for(col = 0; col < dim; col++)
        {
            if ( rig == dim -2 && col == dim -2 )
            {
                filter -> data[liv][rig][col] = 8;
            }
            else
                filter -> data[liv][rig][col] = -1;
        }
    }
    /*copia i valori anche negli altri due livelli*/
    for(rig = 0; rig < filter -> h; rig++)
    {
        for(col = 0; col < filter -> w; col++)
        {
            /*ai canali successivi al primo (lo 0) assegno i valori presenti nel primo. Ovvero, copio i valori
            del primo livello nel secondo e nel terzo, cella per cella, scorrendo i canali. Poi passo all'elemento
            successivo del primo livello (col++) e ripeto*/
            for(liv = 1; liv < filter -> k; liv++)
            {
                filter -> data[liv][rig][col] = filter -> data[0][rig][col];
            }
        }
    }
    compute_stats(filter);

    return filter;
}

/*
CREATE_EMBOSS_FILTER
Autore: Federico
Descrizione: Crea un filtro per aggiungere profondità */
ip_mat * create_emboss_filter()
{
    unsigned int liv, rig, col, dim = 3;
    ip_mat *filter = ip_mat_create(dim, dim, dim, 0);

    if(filter == NULL)
    {
        printf("Errore! Impossibile allocare Ipmat\n");
        exit(1);
    }

    filter -> data[0][0][0] = -2;
    filter -> data[0][0][1] = -1;
    filter -> data[0][0][2] = 0;

    filter -> data[0][1][0] = -1;
    filter -> data[0][1][1] = 1;
    filter -> data[0][1][2] = 1;

    filter -> data[0][2][0] = 0;
    filter -> data[0][2][1] = 1;
    filter -> data[0][2][2] = 2;

    /*copia i valori anche negli altri due livelli*/
    for(rig = 0; rig < filter -> h; rig++)
    {
        for(col = 0; col < filter -> w; col++)
        {
            /*ai canali successivi al primo (lo 0) assegno i valori presenti nel primo. Ovvero, copio i valori
            del primo livello nel secondo e nel terzo, cella per cella, scorrendo i canali. Poi passo all'elemento
            successivo del primo livello (col++) e ripeto*/
            for(liv = 1; liv < filter -> k; liv++)
            {
                filter -> data[liv][rig][col] = filter -> data[0][rig][col];
            }
        }
    }

    compute_stats(filter);
    return filter;
}

/*
CREATE_AVERAGE_FILTER
Autore: Federico
Descrizione:Crea un filtro medio per la rimozione del rumore
*/
ip_mat * create_average_filter(unsigned int h, unsigned int w, unsigned int k)
{
    unsigned int liv, rig, col;
    float c = 1.0/(w*h);
    ip_mat *filter;

    if(h != w)
    {
        printf("Errore! Dimensioni righe e colonne errate!\n");
        exit(1);
    }

    filter = ip_mat_create(h, w, k, 0);

    for(liv = 0; liv < k; liv++)
    {
        for(rig = 0; rig < h; rig++)
        {
            for(col = 0; col < w; col++)
            {
                filter -> data[liv][rig][col] = c;
            }
        }
    }

    compute_stats(filter);
    return filter;
}


/*CREATE_GAUSSIAN_FILTER
Autore: Federico
Descrizione: Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma)
{
    unsigned int liv, rig, col;
    ip_mat *filter;
    filter = NULL;
    filter = ip_mat_create(h, w, k, 0);
    /*per trovare il centro del filtro, se le dimensioni sono uguali, vale da origine sia per righe che per le colonne*/
    /*int centro = w / 2;*/
    /*per la calcolare la somma dei valori del kernel (filter -> data), necessaria per la normalizzazione*/
    double somma[3]; /* abbiamo messo 3 perchè sulle specifiche del progetto è scritto che lavoreremo sempre con matrici data contenenti solo 3 canali per rappresentare il range di colori R G B */
    double x, y;

    for (liv = 0; liv < 3; liv++)
        somma[liv] = 0.0;

    /*se le dimensioni passate non identificano una matrice quadrata*/
    if( h != w || filter == NULL )
    {
        printf("Errore! Dimensioni righe e colonne errate!\n");
        exit(1);
    }



    for(liv = 0; liv < k; liv++)
    {
        for(rig = 0; rig < h; rig++)
        {
            for(col = 0; col < w; col++)
            {
                /*traduzione della formula della distribuzione gaussiana
                filter -> data[liv][rig][col] = exp( -0.5 * ( pow( (rig - centro) / sigma, 2.0 )
                                                + pow( (col - centro) / sigma, 2.0 )) )
                                                    / (2 * M_PI * sigma * sigma);
                filter -> data[liv][rig][col] = ;
                somma[liv] += filter -> data[liv][rig][col];
                */
                x = rig - ( h - 1 ) / 2.0;
                y = col - ( w - 1 ) / 2.0;
                filter->data[liv][rig][col] = exp( (pow(x,2)+pow(y,2)) / (2 * pow(sigma,2) * (-1)) );
                somma[liv] += filter->data[liv][rig][col];
            }
        }
    }

    /*normalizzazione del filtro*/
    for(liv = 0; liv < k; liv++)
    {
        for(rig = 0; rig < h; rig++)
        {
            for(col = 0; col < w; col++)
            {
                filter -> data[liv][rig][col] /= somma[liv];
            }
        }
    }
    return filter;
}

/*
CLAMP
Autore: Tom
Descrizione; controlla che i valori dei tre canali di ogni pixel siano compresi tra 0 e 255 (estremi compresi)
*/
void clamp(ip_mat * t, float low, float high) {
	unsigned int liv, rig, col;

    if(t == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

	for(liv = 0; liv < t->k; liv++) {
    	for(rig = 0; rig < t->h; rig++) {
    		for(col = 0; col < t->w; col++) {
    		    if ( t->data[liv][rig][col] > high ) {
        			t->data[liv][rig][col] = high;
        		}
        		else if ( t->data[liv][rig][col] < low ) {
        			t->data[liv][rig][col] = low;
        		}
			}
		}
	}
}


/*
RESCALE
Autori: Tutto il gruppo in chiamata
Descrizione: Scaliamo la matrice secondo la formula (valore-min)/(max - min). Max e min li prendo da stats e variano da canale a canale.
	     Poi moltiplichiamo i valori della matrice per new_max in modo tale da avere tutti i valori della matrice tra 0 e new_max.
*/
void rescale(ip_mat * t, float new_max) {
	unsigned int liv, rig, col;

    if(t == NULL)
    {
        printf("Errore! IPMAT non valida\n");
        exit(1);
    }

	for(liv = 0; liv < t->k; liv++) {
        for(rig = 0; rig < t->h; rig++) {
        	for(col = 0; col < t->w; col++) {
				/* Scalo la matrice secondo la formula (valore-min)/(max - min). */
				t->data[liv][rig][col] = ( t->data[liv][rig][col] - t->stat[liv].min )  /  ( t->stat[liv].max - t->stat[liv].min );
				/* Moltiplico i valori della matrice per new_max */
				t->data[liv][rig][col] = t->data[liv][rig][col] * new_max;
			}
		}
	}
}

/*
IP_MAT_CONVOLVE
Autore: Tom
Descrizione: Viene fatta la convoluzione della ip_mat *a (applicando il filtro ip_mat *f) su una nuova ip_mat della stessa dimensione
	     di a. I parametri in ingresso non vengono modificati.
*/
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f) {
	ip_mat *x,*y;
	int pad_h, pad_w;
	float value; /* value è la variabile dove mi salvo il valore calcolatoblend dall'incrocio dei valori di a su f */
	unsigned int liv, rig, col, rig1, col1, i, j; /* rig1 e col1 li utilizzerò per muovermi su a->data per ottenere i valori da calcolare sulla matrice filtro */

    if(a == NULL || f == NULL)
    {
        printf("Errore! IPMAT non valida/e\n");
        exit(1);
    }

	/* imposto il valore del pad_h e pad_w
	BISOGNA CONTROLLARE CHE LA MATRICE FILTRO SIA QUADRATA??? BOOOH */
	pad_h = ( f->h - 1 ) / 2;
	pad_w = ( f->w - 1 ) / 2;

	/* creo la nuova ip_mat x su cui scriverò tutte le modifiche */
	y = ip_mat_create(a->h, a->w, a->k, 0);

	/* applico la funzione padding ad x (in realtà non serve perchè con la funzione precedente ho già inizializzato tutta la matrice a 0.0
	   però la funzione a quanto pare bisogna utilizzarla comunque quindi la richiamo lo stesso :) */
	x = ip_mat_padding(y, pad_h, pad_w);
    ip_mat_free(y);
	/* Salvo dentro x->data il risultato del filtro sulla matrice a */
	for(liv = 0; liv < x->k; liv++) {
        	for(rig = pad_h; rig < x->h - pad_h; rig++) {
            		for(col = pad_w; col < x->w - pad_w; col++) {
        				/* ora calcolo il risultato da inserire */
        				value = 0.0; /* resetto value per non avere valori "sporchi" */
                        i = 0;
                        for ( rig1 = rig-pad_h; i < f->h; rig1++, i++ ) {
                            j = 0;
        					for ( col1 = col-pad_w; j < f->w; col1++, j++) {
        						value += ( a->data[liv][rig1][col1] * f->data[liv][i][j] );
        					}
        				}
                x->data[liv][rig][col] = value ;
			}
		}
	}

	return x;
}

/* ####################
   ### FINE PARTE 3 ###
   #################### */

/* ************** FUNZIONI GIÀ IMPLEMENTATE ***********************************
*************************************************************************** */
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

    compute_stats(out);

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
    if(i<a->h && j<a->w &&k<a->k){
        return a->data[k][i][j];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[k][i][j]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}
/* ************** FUNZIONI GIÀ IMPLEMENTATE ***********************************
*************************************************************************** */
