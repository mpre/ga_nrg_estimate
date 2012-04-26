// Include vari
#include <iostream>
#include <cmath>
#include <fstream>

// Inclusione tipo GA
#include <ga/GASimpleGA.h> 
#include <ga/GA1DBinStrGenome.h>
#include <ga/GABin2DecGenome.h>
#include <ga/GASelector.h>

/* Numero di fenotipi del genoma */
const int PHENOTYPE_N   =10;
/* Valori minimi e massimi che i fenotipi possono assumere*/
const int PHENOTYPE_MIN =-10;
const int PHENOTYPE_MAX = 10;

/* Dimensione del torneo di selezione all'interno
   della popolazione */
const int TOURNAMENT_DIM =6;

/* Percentuale di errore oltre la quale viene 
   segnalato un errore (e la linea relativa) 
   all'interno della procedura finale di stampa*/
const int ERR_PERC =35;

/* Numero di giorni precedenti a quello da calcolare
   da prendere in considerazione durante la procedura di 
   fitness */
const int BACKP_NUM =4;

/* Numero di righe e di colonne da leggere nel file */
int nrows = -1;
int nvars = -1;

/* Numero di fitness evaluations
   utilizzato per terminare la computazione */
int ftn_eval =0;

/* Valori letti */
float ** values;

/* Sigma della fitness sharing */
float sh_sg =0.5;

/* Inizializzazione valori tramite lettura da file */
void init_values();

/* Funzione di fitness */
float ftn(GAGenome&);

/* Fitness scaling */
float distance(const GAGenome&, const GAGenome&);

/* Funzioni di utilità (stampa valori etc etc) */
void print_poly(const GAGenome& g);
void ftn_stampata(GAGenome& g);

/* --- MAIN --- */
int main(int argc, char** argv) {
    init_values();
    
    GARandomSeed();
    GABin2DecPhenotype ph;
    int i;
    for (i=1; i<PHENOTYPE_N; ++i)
        ph.add(16, PHENOTYPE_MIN, PHENOTYPE_MAX);
    GABin2DecGenome genome(ph, ftn);
    GASimpleGA ga(genome);

    GASharing gs(distance, sh_sg);
    
    ga.scaling(gs);

    ga.parameters("GAconfiguration.conf");
    ga.minimize();
    ga.selector(GATournamentSelector(TOURNAMENT_DIM));
    i =0;
    while(ftn_eval < 50000-ga.populationSize()) {
        // Evoluzione
        ++ga;
        // Print miglior individuo
        std::cout << "-- Evoluzione #" << (i+1) << std::endl;
        std::cout << "-- Numero di valutazioni di fitness effettuate: " << ftn_eval << std::endl;
        std::cout << "-- Miglior individuo:  " << ga.statistics().bestIndividual() << std::endl;
        std::cout << "Fitness: " << ftn((const_cast<GAGenome&>(ga.statistics().bestIndividual()))) << std::endl << std::endl;
        ++i;
    }

    ftn_stampata(const_cast<GAGenome&>(ga.statistics().bestIndividual()));
    
    std::cout << "Miglior individuo: ";
    print_poly(ga.statistics().bestIndividual());
    std::cout << "Miglior fitness  : " << std::endl << ftn(const_cast<GAGenome&>(ga.statistics().bestIndividual())) << std::endl;

    for (int i=0; i<nrows; i++)
      delete values[i];
    
}

void print_poly(const GAGenome& g)
{
  GABin2DecGenome& gb = (GABin2DecGenome&) g;
  int ncoeff = gb.nPhenotypes();
  for (int j = 0; j < ncoeff; j++) {
    std::cout << gb.phenotype(j) <<"\t";
  }
  std::cout << std::endl;
}

void init_values() {
    std::ifstream f;
    f.open("energia.txt");
    f >> nrows;
    f >> nvars;
    values = new float*[nrows];
    for(int i=0; i<nrows; ++i) {
      values[i] = new float[nvars];
      for(int j=0; j<nvars; ++j) {
        f >> values[i][j];
      }
    }
    f.close();
}

float ftn(GAGenome& g) {
  ftn_eval++;
  GABin2DecGenome& gb = (GABin2DecGenome&) g;
  int ncoeff = gb.nPhenotypes();
  float score = 0.0f;
  for (int i = 0; i < nrows; i++) {
    float val = 0;
    for (int j = 7-BACKP_NUM; j < ncoeff; j++) {
      float tmp =0;

      if(j == 8) {
        tmp += 1400 * (gb.phenotype(8)*values[i][8] - gb.phenotype(7)*values[i][7] );
      }
      if(j>=7-BACKP_NUM && j<7) {
        tmp += (gb.phenotype(j)*values[i][j]/(std::log(7-j)+1) - gb.phenotype(j-1)*values[i][j-1]/(std::log(7-j+1)+1));
      }

      val += tmp;
    }
    score += std::pow(val - values[i][9],2);
  }
  score = std::sqrt(score/nrows);
  return score;

}

void ftn_stampata(GAGenome& g) {
  /* Questa funzione serve solamente a stampare alla 
   * fine della procedura di calcolo i valori previsti per
   * ogni caso e il loro scostamento dal reale.
   */
  int err_count=0;
  GABin2DecGenome& gb = (GABin2DecGenome&) g;
  int ncoeff = gb.nPhenotypes();
  float score = 0.0f;
  for (int i = 0; i < nrows; i++) {
    float val = 0;
    for (int j = 0; j < ncoeff; j++) {
      float tmp =0;

      if(j == 8) {
        tmp += 1400 * (gb.phenotype(8)*values[i][8] - gb.phenotype(7)*values[i][7] );
      }
      if(j>=7-BACKP_NUM && j<7) {
        tmp += (gb.phenotype(j)*values[i][j]/(std::log(7-j)+1) - gb.phenotype(j-1)*values[i][j-1]/(std::log(7-j+1)+1));
      }

      val += tmp;
    }
    score += std::pow(val - values[i][9],2);
    std::cout <<i+1<<"->\t" << "Calcolato: " << val << "\tReale: " << values[i][9];
    std::cout << "\tDifferenza: " << val-values[i][9]<< "\tErrore %: " << 100*((val - values[i][9])/val) <<"%" ;
    
    if(abs(100*((val - values[i][9])/val)) > ERR_PERC) {
      std::cout << "\tFile line #"<<i+3<<" !";
      err_count++;
    }
    std::cout << std::endl;
  }

  score = std::sqrt(score/nrows);

  std::cout << std::endl << std::endl;
  std::cout << "Numero di valutazioni di fitness effettuate: " << ftn_eval << std::endl;
  std::cout << "Risultati con percentuale di errore maggiore di ";
  std::cout << ERR_PERC << ": " << err_count << " / " << nrows << std::endl;

}


float distance(const GAGenome& g1, const GAGenome& g2) {
    GABin2DecGenome& gb1 = (GABin2DecGenome&) g1;
    GABin2DecGenome& gb2 = (GABin2DecGenome&) g2;
    
    float dist =0.0f;
    for(int i=0; i<gb1.length(); i++) {
      if(abs(gb1.gene(i) - gb2.gene(i)) < 3.f) {
          dist++;
        }
    }
    return dist/gb1.length();
}
