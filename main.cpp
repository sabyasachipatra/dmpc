#include <bits/stdc++.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
// using namespace std;
#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "Graph.h"
#include "Subgraph.h"
#include "SubgraphTree.h"
#include "GraphIsomor.h"
#include "Timer.h"


#define MAX_BUF 256  // Maximum string buffer size
#define MAX_ITERATION_COUNT 30   // Maximum number of iterations

struct edgelist
{
    int n;
    int level;
    struct edgelist *next;
};


void parse_cmdline(int argc, char **argv);
void check_inputs();
void initialize();
void prepare_graph();
void prepare_testgraph();
void compute_features();
void train_model();
void predict_complex();
void extract_cliques();
void extract_clique(int a);
void read_benchmark();
void extractFeature(Graph* gc,double feature[],double norm_feature[]);
float correlationCoefficient(int X[], int Y[], int n);
int catagorize(double v);
// Internally Defined Routines
int Singular_Value_Decomposition(double* A, int nrows, int ncols, double* U, double* singular_values, double* V, double* dummy_array);
static void Householders_Reduction_to_Bidiagonal_Form(double* A, int nrows, int ncols, double* U, double* V, double* diagonal, double* superdiagonal );
static int  Givens_Reduction_to_Diagonal_Form( int nrows, int ncols, double* U, double* V, double* diagonal, double* superdiagonal );
static void Sort_by_Decreasing_Singular_Values(int nrows, int ncols, double* singular_value, double* U, double* V);
void floydWarshall(int matrix[30][30], int nV);

void printedgelist(struct edgelist *r);
struct edgelist* insertedge(struct edgelist *r, int v);
struct edgelist* insertsortedge(struct edgelist *r, int v);
struct edgelist* deleteedge(struct edgelist *r, int v);
void printedgelistwitha(struct edgelist *r, int a);
void insertneighbor(std::vector<int> *nei,int a, struct edgelist *temp);
void setlevel(struct edgelist *r, int l);
int duplicate(struct edgelist *r);
void print_complex();
void output_complex();
int compute_lengths(struct edgelist *r);
int subcomplex(struct edgelist *ea);

// Variable declarations
int num_exchanges;
int num_tries;
char graph_file[MAX_BUF];
char testgraph_file[MAX_BUF];
char map_file[MAX_BUF];
char testmap_file[MAX_BUF];
char method[MAX_BUF];
char benchmark_file[MAX_BUF];
Graph *g;
Graph *testg;
std::vector<Subgraph*> sgv;
std::vector<Subgraph*> sgv_random;
SubgraphTree isomorSG;
std::vector<Graph*> graphv;

std::vector<edgelist*> complexes;
std::vector<edgelist*> out3complexes;
std::vector<edgelist*> outcomplexes;
std::vector<int> lengths;
std::vector<edgelist*> newcomplexes;
std::vector<edgelist*> addedcomplexes;



int main(int argc, char **argv) {
    std::cout << "This program is developed by Sabyasachi Patra, ";
    std::cout << "IIIT Bhubaneswar, India" << std::endl;
    initialize();
    parse_cmdline(argc, argv);
    check_inputs();
    prepare_graph();
    prepare_testgraph();
    compute_features();
    train_model();
    predict_complex();
    print_complex();
    output_complex();
    GraphIsomor::finishNauty();
    return 0;
}

// Initialize everything
void initialize() {
    num_exchanges = 2;
    num_tries = 7;
    snprintf(method, strlen("none")+1, "%s\n", "none");
    snprintf(graph_file, strlen("none")+1, "%s\n", "none");
    snprintf(benchmark_file, strlen("none")+1, "%s\n", "none");
    FILE *fp;
    fp = fopen("Output.txt", "w");
    fprintf(fp, "OUTPUT");
    fprintf(fp, "\n");
    fclose(fp);
    fp = fopen("Error.txt", "w");
    fprintf(fp, "ERROR");
    fprintf(fp, "\n");
    fclose(fp);
    GraphIsomor::initNauty(10);
}

// Parse all command line arguments
void parse_cmdline(int argc, char **argv) {
    for (int i=1; i < argc; i++) {
        // Graph file
        if (!strcmp("-g", argv[i]) || !strcmp("--graph", argv[i])) {
            snprintf(graph_file, strlen(argv[i+1])+1, "%s\n", argv[i+1]);
            // cout << "graph_file = " << graph_file << endl;
            i++;
        }
        else if (!strcmp("-t", argv[i]) || !strcmp("--testgraph", argv[i])) {
            snprintf(testgraph_file, strlen(argv[i+1])+1, "%s\n", argv[i+1]);
            // cout << "testgraph_file = " << testgraph_file << endl;
            i++;
        }
        else if (!strcmp("-m", argv[i]) || !strcmp("--method", argv[i])) {
            snprintf(method, strlen(argv[i+1])+1, "%s\n", argv[i+1]);
            // cout << "method = " << method << endl;
            if (strcmp(method, "benchmark") == 0) {
                    strcpy(benchmark_file, argv[i+2]);
                    i++;
            }
            i++;
        }
    }
}

void check_inputs() {
    if (strcmp(method, "benchmark") != 0) {
        std::cout << "invalid method" << std::endl;
        exit(1);
    }

    if (strcmp(benchmark_file, "none") == 0) {
        //std::cout << "no input graph file" << std::endl;
        exit(1);
    }

    if (strcmp(graph_file, "none") == 0) {
        //std::cout << "no input graph file" << std::endl;
        exit(1);
    }
}

// Prepare the real graph for computation
void prepare_graph() {
    g = new Graph();
    // Read the graph file
    strcat(map_file, "map_");
    strcat(map_file, graph_file);
    Graph::readGraphFile(g, graph_file, map_file);
    // sort and create neighbours array
    g->sortNeighbours();
    g->makeNeighboursArray();
    printf("graph file: %s\n", graph_file);
    printf("%d nodes, %d edges\n", g->numberNodes(), g->numberEdges());
}

void prepare_testgraph() {
    testg = new Graph();
    // Read the graph file
    strcat(testmap_file, "map_");
    strcat(testmap_file, testgraph_file);
    Graph::readGraphFile(testg, testgraph_file, testmap_file);
    // sort and create neighbours array
    testg->sortNeighbours();
    testg->makeNeighboursArray();
    printf("testgraph file: %s\n", testgraph_file);
    printf("%d nodes, %d edges\n", testg->numberNodes(), testg->numberEdges());
}


// Compute features from benchmark complexes
void compute_features() {
  int i,j, cntr=0;
  Timer::start(0);
  read_benchmark();
  Timer::stop(0);
  printf("Creation time: %.2f\n", Timer::elapsed(0));
  std::vector<Subgraph *>::iterator ii;
  std::cout << "sgv.size() = " << sgv.size() << std::endl;

  std::ofstream outfile;
  std::ofstream normfile;
  outfile.open("features.csv");
  outfile << "NodeSize,GraphDensity,meanDegree,varDegree,medianDegree,maxDegree,";
  outfile << "meanDegreeCorrelation,varDegreeCorrelation,maxDegreeCorrelation,";
  outfile << "meanClusteringCoeff,varClusteringCoeff,maxClusteringCoeff,";
  outfile << "meanTopologicCoeff,varTopologicCoeff,maxTopologicCoeff,";
  outfile << "eigenValue1,eigenValue2,eigenValue3,";
  outfile << "aveLength,maxLength,aveWeight,maxWeight,level\n";

  normfile.open("normfeatures.csv");
  normfile << "NodeSize,GraphDensity,meanDegree,varDegree,medianDegree,maxDegree,";
  normfile << "meanDegreeCorrelation,varDegreeCorrelation,maxDegreeCorrelation,";
  normfile << "meanClusteringCoeff,varClusteringCoeff,maxClusteringCoeff,";
  normfile << "meanTopologicCoeff,varTopologicCoeff,maxTopologicCoeff,";
  normfile << "eigenValue1,eigenValue2,eigenValue3,";
  normfile << "aveLength,maxLength,aveWeight,maxWeight,level\n";

  std::ofstream randoutfile;
  std::ofstream randnormfile;
  randoutfile.open("randfeatures.csv");
  randoutfile << "NodeSize,GraphDensity,meanDegree,varDegree,medianDegree,maxDegree,";
  randoutfile << "meanDegreeCorrelation,varDegreeCorrelation,maxDegreeCorrelation,";
  randoutfile << "meanClusteringCoeff,varClusteringCoeff,maxClusteringCoeff,";
  randoutfile << "meanTopologicCoeff,varTopologicCoeff,maxTopologicCoeff,";
  randoutfile << "eigenValue1,eigenValue2,eigenValue3,";
  randoutfile << "aveLength,maxLength,aveWeight,maxWeight,level\n";

  randnormfile.open("randnormfeatures.csv");
  randnormfile << "NodeSize,GraphDensity,meanDegree,varDegree,medianDegree,maxDegree,";
  randnormfile << "meanDegreeCorrelation,varDegreeCorrelation,maxDegreeCorrelation,";
  randnormfile << "meanClusteringCoeff,varClusteringCoeff,maxClusteringCoeff,";
  randnormfile << "meanTopologicCoeff,varTopologicCoeff,maxTopologicCoeff,";
  randnormfile << "eigenValue1,eigenValue2,eigenValue3,";
  randnormfile << "aveLength,maxLength,aveWeight,maxWeight,level\n";


  std::ofstream trainfile;
  trainfile.open("trainfeatures.csv");
  trainfile << "NodeSize,GraphDensity,meanDegree,varDegree,medianDegree,maxDegree,";
  trainfile << "meanDegreeCorrelation,varDegreeCorrelation,maxDegreeCorrelation,";
  trainfile << "meanClusteringCoeff,varClusteringCoeff,maxClusteringCoeff,";
  trainfile << "meanTopologicCoeff,varTopologicCoeff,maxTopologicCoeff,";
  trainfile << "eigenValue1,eigenValue2,eigenValue3,";
  trainfile << "aveLength,maxLength,aveWeight,maxWeight,level\n";


  for(ii=sgv.begin(); ii!=sgv.end(); ii++) {

	//std::cout << (*ii)->canstr << std::endl;
    Graph *gc = new Graph();
	int sz=sqrt(strlen((*ii)->canstr));
    Graph::strToGraph(gc, (*ii)->canstr, sz);
    gc->sortNeighbours();
    gc->makeNeighboursArray();
    double feature[22];
    double norm_feature[22];
    if(gc->numberNodes() < 3) continue;
    if(gc->numberEdges() == 0) continue;
    extractFeature(gc,feature, norm_feature);
    for(i=0;i<22;i++) {
        //std::cout<< std::setw(12) << feature[i] << "  ";
        outfile << feature[i] << ",";
    }
    //td::cout << std::endl;
    outfile << "yes\n";

    for(i=0;i<22;i++) {
        //std::cout<< std::setw(12) << norm_feature[i] << "  ";
        normfile << norm_feature[i] << ",";
    }
    //std::cout << std::endl;
    normfile << "yes\n";

    gc->deleteGraph();
    cntr++;
  }
  outfile.close();
  normfile.close();
  int rc=0, rnc=0;
  for(ii=sgv_random.begin(); ii!=sgv_random.end(); ii++) {

	//std::cout << (*ii)->canstr << std::endl;
    Graph *gc = new Graph();
	int sz=sqrt(strlen((*ii)->canstr));
    Graph::strToGraph(gc, (*ii)->canstr, sz);
    gc->sortNeighbours();
    gc->makeNeighboursArray();
    double randfeature[22];
    double norm_randfeature[22];
    if(gc->numberNodes() < 3) continue;
    extractFeature(gc,randfeature,norm_randfeature);
    for(i=0;i<22;i++) {
        //std::cout<< std::setw(12) << randfeature[i] << "  ";
        randoutfile << randfeature[i] << ",";
    }
    //std::cout << std::endl;
    randoutfile << "no\n";

    for(i=0;i<22;i++) {
        //std::cout<< std::setw(12) << norm_randfeature[i] << "  ";
        randnormfile << norm_randfeature[i] << ",";
    }
    //std::cout << std::endl;
    randnormfile << "no\n";

    gc->deleteGraph();
    rc++;
    rnc++;
  }
  printf("rc=%d    rnc=%d\n",rc,rnc);
  randoutfile.close();
  randnormfile.close();

  std::string line1;
  std::string line2;
  std::ifstream f1 ("normfeatures.csv");
  std::ifstream f2 ("randnormfeatures.csv");
  if (!f1) {
    std::cout << "invalid file" << std::endl;
    exit(1);
  }
  if (!f2) {
    std::cout << "invalid file" << std::endl;
    exit(1);
  }
  getline (f1,line1);
  getline (f2,line2);
  cntr=0;
  while ( getline (f1,line1)) {
     trainfile << line1 << "\n";
     if (getline (f2,line2)) trainfile << line2 << "\n";
     if (getline (f2,line2)) trainfile << line2 << "\n";
     if (getline (f2,line2)) trainfile << line2 << "\n";
     if (getline (f2,line2)) trainfile << line2 << "\n";
     //if (getline (f2,line2)) trainfile << line2 << "\n";
  }
  trainfile.close();
  f1.close();
  f2.close();

}

void train_model()
{
    printf("training start\n");
    system("py train.py >> train_out.txt");
}

void output_complex()
{
  std::string line;
  std::ifstream mfile(map_file);
  std::cout << map_file << "\n";
  if (!mfile) {
    std::cout << "invalid file" << std::endl;
    exit(1);
  }
  std::map<int, std::string> maps;
  while ( getline (mfile,line) ) {
     //std::cout << line << '\n';
     // Vector of string to save tokens
     std::vector <std::string> tokens;
     // stringstream class check1
     std::stringstream check1(line);
     std::string intermediate;
     // Tokenizing w.r.t. space ' '
     while(std::getline(check1, intermediate, '\t'))
     {
         tokens.push_back(intermediate);
     }
     //std::cout << "tokens.size= " << tokens.size() << "\n";
     maps.insert(std::pair<int, std::string>(stoi(tokens[1]), tokens[0]));
  }
//  std::map<int, std::string>::iterator itr;
//  for (itr=maps.begin(); itr != maps.end(); ++itr) {
//      std::cout << itr->first << '\t' << itr->second << "\n";
//  }

    std::vector<edgelist*>::iterator ii;
    printf(".......................\n\n");
    std::ofstream f1 ("outputcomplexes.txt");
        if (!f1) {
            std::cout << "invalid file" << std::endl;
            exit(1);
        }
    for(ii=outcomplexes.begin(); ii!=outcomplexes.end();ii++) {
        printedgelist(*ii);
        struct edgelist *temp=(*ii);
        while(temp != NULL) {
            std::cout << maps[temp->n] << "  ";
            if(temp->next != NULL) f1 << maps[temp->n] << ",";
            else f1 << maps[temp->n];
            temp=temp->next;
        }
        f1 << "\n";
        printf("\n");
    }
    f1.close();
    printf(".......................\n\n");
}

void print_complex()
{
    std::vector<edgelist*>::iterator ii;
    std::vector<int>::iterator jj;
    //std::vector<edgelist*> complexes;
    std::vector<int> lengths;
    for(ii=complexes.begin(); ii!=complexes.end();ii++) {
        int l=compute_lengths(*ii);
        //printf("%d ",l);
        lengths.push_back(l);
    }
//    for(jj=lengths.begin(); jj!=lengths.end();jj++) {
//        printf("%d ",*jj);
//    }
//    printf("\n");
    while(lengths.size()>0) {
        int i=0,j=0;
        int maxl=0;
        for(jj=lengths.begin(); jj!=lengths.end();j++,jj++) {
            if(*jj > maxl) {
                maxl=*jj;
                i=j;
            }
        }
        //printf("l=%d  max=%d\n",lengths.size(),maxl);
        for(ii=complexes.begin(),j=0; ii!=complexes.end();ii++,j++) {
            if(j==i) break;
        }
        if (subcomplex(*ii)==0) {
            //printf("new complex\n");
            outcomplexes.push_back(*ii);
        }
        //else printf("sub complex\n");
        lengths.erase(lengths.begin()+i);
        complexes.erase(complexes.begin()+i);
    }
    //for(ii=out3complexes.begin(); ii!=out3complexes.end();ii++) {
    //    printf("size 3\n");
    //    outcomplexes.push_back(*ii);
    //}
    //for(ii=outcomplexes.begin(); ii!=outcomplexes.end();ii++) {
    //    printedgelist(*ii);
    //}
}

int subcomplex(struct edgelist *ea)
{
    std::vector<edgelist*>::iterator ii;
    for(ii=outcomplexes.begin(); ii!=outcomplexes.end();ii++) {
        struct edgelist *temp=(*ii);
        while(ea != NULL) {
            while(temp != NULL) {
                if(ea->n == temp->n) break;
                temp=temp->next;
            }
            if (temp == NULL) break;
            ea=ea->next;
        }
        if (ea == NULL) {
            return 1;
        }
    }
    return 0;
}

int compute_lengths(struct edgelist *ea)
{
    int l=0;
    while(ea != NULL) {
        ea=ea->next;
        l++;
    }
    return l;
}

void predict_complex()
{
    extract_cliques();

    std::vector<edgelist*>::iterator ii;
    int i,j,k,start=0, lvl=0;
    printf(".......................\n\n");
    for(ii=complexes.begin(); ii!=complexes.end();ii++) {
        newcomplexes.push_back(*ii);
        //printedgelist(*ii);
    }
    printf(".......................\n\n");
    while(1) {
        for(ii=newcomplexes.begin(); ii!=newcomplexes.end(); ii++) {
            struct edgelist *temp=(*ii);
            std::vector<int> nei;
            int gsize=0;
            while(temp != NULL) {
                int *na = testg->getNeighboursArray(temp->n);
                for(j=0;j<testg->numberNeighbours(temp->n);j++) {
                    insertneighbor(&nei,na[j],*ii);
                }
                temp=temp->next;
                gsize++;
            }
            std::vector<int>::iterator jj;
//            for(jj=nei.begin(); jj!=nei.end(); jj++) {
//                printf("%d  ", *jj);
//            }
//            printf("\n\n");

            temp=(*ii);
            std::map<int, int> nmap;
            j=0;
            while(temp != NULL) {
                nmap.insert(std::pair<int, int>(j,temp->n));
                temp=temp->next;
                j++;
            }

//            while(nei.size()>0) {
//                //printf("nei.size=%d\n",nei.size());
                int newsize=gsize;
                int flag=0, l=0;
                std::vector<int> nnei;
                for(jj=nei.begin(); jj!=nei.end(); jj++) {
                    //printf("jj=%d\n", *jj);
                    nmap.insert(std::pair<int, int>(newsize,*jj));
                    newsize++;
                    //printf("newsize=%d\n",newsize);
                    Graph *gc = new Graph();
                    gc->createGraph(newsize);
                    for (j=0; j < newsize; j++) {
                        for (k=0; k < newsize; k++) {
                            if (testg->hasEdge(nmap[j], nmap[k])) {
                                gc->addEdge(j,k);
                            }
                        }
                    }
                    gc->sortNeighbours();
                    gc->makeNeighboursArray();
                    int n=gc->numberNodes();
                    int e=gc->numberEdges();
                    double dens=2.0*e/(n*(n+1));
                    //printf("e=%d  dens(%d)=%lf\n",e,n,dens);
                    if (dens<0.6) {
                        nmap.erase(newsize-1);
                        newsize--;
                        if(flag==1) {
                            //nnei.push_back(*jj);
                            struct edgelist *ea=NULL;
                            for(j=0;j<newsize;j++) {
                                //printf("j=%d nmap[j]=%d  ",j,nmap[j]);
                                ea=insertsortedge(ea,nmap[j]);
                            }
                            //setlevel(ea,lvl+1);
                            //printf("\nComplex added\n");
                            //printedgelist(ea);

                            addedcomplexes.push_back(ea);
                            for(j=gsize;j<newsize;j++) {
                                nmap.erase(j);
                            }
                            flag=0;
                            newsize=gsize;
                        }
                        continue;
                    }
                    else {
                        flag=1;
                    }

                }
//                nei.clear();
//                for(jj=nnei.begin(); jj!=nnei.end(); jj++) {
//                    nei.push_back(*jj);
//                }
//                //printf("nei.size=%d\n",nei.size());
//            }
        }
        newcomplexes.clear();
        printf("addedcomplexes.size=%d\n",addedcomplexes.size());
        std::ofstream testfile;
            testfile.open("testfeatures.csv");
            testfile << "NodeSize,GraphDensity,meanDegree,varDegree,medianDegree,maxDegree,";
            testfile << "meanDegreeCorrelation,varDegreeCorrelation,maxDegreeCorrelation,";
            testfile << "meanClusteringCoeff,varClusteringCoeff,maxClusteringCoeff,";
            testfile << "meanTopologicCoeff,varTopologicCoeff,maxTopologicCoeff,";
            testfile << "eigenValue1,eigenValue2,eigenValue3,";
            testfile << "aveLength,maxLength,aveWeight,maxWeight,level\n";

        for(ii=addedcomplexes.begin(); ii!=addedcomplexes.end();ii++) {
            struct edgelist *temp=(*ii);
            std::map<int, int> nmap;
            j=0;
            while(temp != NULL) {
                nmap.insert(std::pair<int, int>(j,temp->n));
                temp=temp->next;
                j++;
            }
            int newsize=j;
            Graph *gc = new Graph();
            gc->createGraph(newsize);
            for (j=0; j < newsize; j++) {
                for (k=0; k < newsize; k++) {
                    if (testg->hasEdge(nmap[j], nmap[k])) {
                        gc->addEdge(j,k);
                    }
                }
            }
            gc->sortNeighbours();
            gc->makeNeighboursArray();
            int n=gc->numberNodes();
            int e=gc->numberEdges();
            double feature[22];
            double norm_feature[22];
            if(gc->numberNodes() < 3) continue;
            if(gc->numberEdges() == 0) continue;
            extractFeature(gc,feature, norm_feature);
            for(i=0;i<22;i++) {
                //std::cout<< std::setw(12) << norm_feature[i] << "  ";
                testfile << norm_feature[i] << ",";
            }
            //std::cout << std::endl;
            testfile << "predict\n";
        }
        testfile.close();
        printf("test start\n");
        system("py test.py >> test_out.txt");
        std::string line1;
        std::ifstream f1 ("out.txt");
        if (!f1) {
            std::cout << "invalid file" << std::endl;
            exit(1);
        }
        for(ii=addedcomplexes.begin(); ii!=addedcomplexes.end();ii++) {
            if(duplicate(*ii)==1) continue;
            getline (f1,line1);
            if(line1.compare("no")==0) continue;
            //std::cout << line1 << endl;
            newcomplexes.push_back(*ii);
            complexes.push_back(*ii);
        }
        f1.close();
        addedcomplexes.clear();
        printf("newcomplexes.size=%d\n",newcomplexes.size());
        if(newcomplexes.size()<=0) break;
    }
//    printf(".......................\n\n");
//    for(ii=complexes.begin(); ii!=complexes.end();ii++) {
//        printedgelist(*ii);
//    }
//    printf(".......................\n\n");
}

void insertneighbor(std::vector<int> *nei,int a, struct edgelist *temp)
{
    while(temp != NULL) {
        if(temp->n==a) {
            //printf("return\n");
            return;
        }
        temp=temp->next;

    }
    std::vector<int>::iterator jj;
    for(jj=(*nei).begin(); jj!=(*nei).end(); jj++) {
        if(*jj==a) return;
    }
    nei->push_back(a);
}
void extract_cliques()
{
    int i,j,k;
    int n=testg->numberNodes();
    int e=testg->numberEdges();
    for(i=1;i<n;i++) {
        printf("i=%d\n",i);
        extract_clique(i);
    }
    std::vector<edgelist*>::iterator ii;

    //for(ii=complexes.begin(); ii!=complexes.end(); ii++) {
    //    printedgelist((*ii));
    //}
}
void extract_clique(int a)
{
    int i,j,k,e=0;
    int n=testg->numberNeighbours(a);
    //printf("n=%d\n",n);
    if(n<2) return;
    e=e+n;
    int *na = testg->getNeighboursArray(a);
    int d[n]; // degree of n neighbors
    int v[n]; // neiigbor
    //int vv[n][n];
    struct edgelist *ea=NULL;
    for (i=0;i<n;i++) ea=insertedge(ea,na[i]);
    struct edgelist *vv[n];
    for (i=0;i<n;i++) vv[i]=NULL;
    for(i=0;i<n;i++) {
        v[i]=na[i];
        int ni=testg->numberNeighbours(na[i]);
        int *nai = testg->getNeighboursArray(na[i]);
        int deg=1; // na[i] is connected to a
        //vv[i][0]=a;
        vv[i]=insertedge(vv[i],a);
        for(j=0;j<ni;j++) {
            if (nai[j]==a) continue;
            for(k=0;k<n;k++) {
                if (nai[j]==na[k]) {
                    //vv[i][deg]=na[k];
                    deg++;
                    vv[i]=insertedge(vv[i],na[k]);
                }
            }
        }
        d[i]=deg;
    }
    for(i=0;i<n;i++) {
        e=e+d[i];
        //printf("d[%d]=%dfrom%d\t",na[i],d[i],g->numberNeighbours(na[i]));
        //printedgelist(vv[i]);
    }
    //printf("\n");
    //printf("e=%d\n",e);
    e=e/2;
    double dens=2.0*e/(n*(n+1));
    //printf("e=%d  dens(%d)=%lf\n",e,n,dens);
    if(n==2 && dens<1) return;

    //all vertices
    //for(i=0;i<n;i++) {
    //    printf("%d  ",v[i]);
    //}
    //printf("\n");

    //sort vertices by the order of their degree
    for(i=0;i<n-1;i++) {
        for(j=n-1;j>i;j--) {
            if(d[j]>d[j-1]) {
                int t=d[j];
                d[j]=d[j-1];
                d[j-1]=t;
                t=v[j];
                v[j]=v[j-1];
                v[j-1]=t;
                struct edgelist *temp;
                temp=vv[j];
                vv[j]=vv[j-1];
                vv[j-1]=temp;
            }
        }
    }
    std::map<int, int> nmap;
    for(i=0;i<n;i++) {
        nmap.insert(std::pair<int, int>(v[i], i));
        //printf("d[%d]=%dfrom%d\t",v[i],d[i],g->numberNeighbours(na[i]));
        //printedgelist(vv[i]);
    }
    //printf("\n");

    while(dens<1 && n>=2) {
        e=e-d[n-1];
        //printedgelist(ea);
        ea=deleteedge(ea,v[n-1]);
        vv[n-1]=deleteedge(vv[n-1],vv[n-1]->n);
        d[n-1]=-1;
        //printedgelist(ea);
        struct edgelist *tedge=vv[n-1];
        while(tedge!=NULL) {
            int index=nmap[tedge->n], nextindex,t;
            vv[index]=deleteedge(vv[index],v[n-1]);
            d[index]--;
            nextindex=index+1;
            while(d[index]<d[nextindex]) {
                t=d[index];
                d[index]=d[nextindex];
                d[nextindex]=t;

                t=v[index];
                v[index]=v[nextindex];
                v[nextindex]=t;

                nmap.erase(v[index]);
                nmap.erase(v[nextindex]);
                nmap.insert(std::pair<int, int>(v[index], index));
                nmap.insert(std::pair<int, int>(v[nextindex], nextindex));

                struct edgelist *temp;
                temp=vv[index];
                vv[index]=vv[nextindex];
                vv[nextindex]=temp;

                nextindex++;
            }

            tedge=tedge->next;
        }
        //printedgelist(vv[n-1]);
        while(vv[n-1]!=NULL) {
            vv[n-1]=deleteedge(vv[n-1],vv[n-1]->n);
        }
        //printedgelist(vv[n-1]);
        n--;
        dens=2.0*e/(n*(n+1));
        //printf("e=%d  dens(%d)=%lf\n",e,n+1,dens);
    }
    if (n>=2) {
        printf("e=%d  dens(%d)=%lf\n",e,n+1,dens);
        //printf("%d ",a);
        //printedgelistwitha(ea,a);
        ea=insertsortedge(ea,a);
        //printedgelist(ea);
        setlevel(ea,0);
        if(duplicate(ea)==0) {
            complexes.push_back(ea);
            if (compute_lengths(ea) == 3) out3complexes.push_back(ea);
        }
    }
}
void setlevel(struct edgelist *r, int l)
{
    while(r!=NULL) {
        r->level=l;
        r=r->next;
    }
}
int duplicate(struct edgelist *ea)
{
    std::vector<edgelist*>::iterator ii;
    int flag=0;
    if (complexes.size()==0) return 0;
    for(ii=complexes.begin(); ii!=complexes.end(); ii++) {
        struct edgelist *r=ea;
        struct edgelist *p=*ii;
        while(r!=NULL && p!= NULL) {
            if(r->n == p->n) {
                r=r->next;
                p=p->next;
            }
            else {
                flag=1;
                break;
            }
        }

        if (flag==0 && r==NULL && p==NULL) {
            //printf("duplicate\n");
            return 1;
        }
        flag=0;
    }
    return 0;
}

struct edgelist* insertedge(struct edgelist *r, int v)
{
    struct edgelist *t;
    struct edgelist *ne;
    ne=(struct edgelist*)malloc(sizeof(struct edgelist));
    ne->n=v;
    ne->next=NULL;
    if (r==NULL) {
        r=ne;
        return r;
    }
    t=r;
    while(t->next!=NULL) t=t->next;
    t->next=ne;
    return r;

}

struct edgelist* insertsortedge(struct edgelist *r, int v)
{
    struct edgelist *t;
    struct edgelist *ne;
    ne=(struct edgelist*)malloc(sizeof(struct edgelist));
    ne->n=v;
    ne->next=NULL;
    if (r==NULL) {
        r=ne;
        return r;
    }
    if (r->n > v) {
        ne->next=r;
        r=ne;
        return r;
    }
    t=r;
    while(t->next!=NULL && t->next->n < v) t=t->next;
    ne->next=t->next;
    t->next=ne;
    return r;
}

struct edgelist* deleteedge(struct edgelist *r, int v)
{
    struct edgelist *t, *t2;

    if (r==NULL) {
        return r;
    }
    if (r->n==v) {
        t=r;
        r=r->next;
        free(t);
        return r;
    }
    int flag=0;
    t=r;
    while(t->next!=NULL) {
        if(t->next->n == v) {
            flag=1;
            break;
        }
        t=t->next;
    }
    if (flag==0)  printf("%d not found\n",v);
    else {
        t2=t->next;
        t->next=t->next->next;
        free(t2);
    }
    return r;

}

void printedgelist(struct edgelist *r)
{
    if (r==NULL) printf("empty list\n");
    while(r!=NULL) {
        printf("%d  ",r->n);
        r=r->next;
    }
    printf("\n");
}

void printedgelistwitha(struct edgelist *r, int a)
{
    int flag=0;
    if (r==NULL) printf("empty list\n");
    while(r!=NULL) {
        if (r->n < a || flag==1) printf("%d  ",r->n);
        else {
            printf("%d  ",a);
            printf("%d  ",r->n);
            flag=1;
        }
        r=r->next;
    }
    if (flag==0) printf("%d  ",a);
    printf("\n");
}

// Populate subgraph vector with subgraphs of 'size' read from file 'subgraphs_file'
void read_benchmark() {
  std::string line;
  char buf[MAX_BUF];
  std::ifstream f (benchmark_file);
  std::cout << benchmark_file << "\n";
  if (!f) {
    std::cout << "invalid file" << std::endl;
    exit(1);
  }

  std::ifstream mfile(map_file);
  std::cout << map_file << "\n";
  if (!mfile) {
    std::cout << "invalid file" << std::endl;
    exit(1);
  }
  std::map<std::string, int> maps;
  while ( getline (mfile,line) ) {
     //std::cout << line << '\n';
     // Vector of string to save tokens
     std::vector <std::string> tokens;
     // stringstream class check1
     std::stringstream check1(line);
     std::string intermediate;
     // Tokenizing w.r.t. space ' '
     while(std::getline(check1, intermediate, '\t'))
     {
         tokens.push_back(intermediate);
     }
     //std::cout << "tokens.size= " << tokens.size() << "\n";
     maps.insert(std::pair<std::string, int>(tokens[0], stoi(tokens[1])));
  }
//  std::map<std::string, int>::iterator itr;
//  for (itr=maps.begin(); itr != maps.end(); ++itr) {
//      std::cout << itr->first << '\t' << itr->second << "\n";
//  }

  while ( std::getline(f,line) ) {
     int i, j, k;
     //std::cout << line << '\n';
     // Vector of string to save tokens
     std::vector <std::string> tokens;
     // stringstream class check1
     std::stringstream check1(line);
     std::string intermediate;
     // Tokenizing w.r.t. space ' '
     while(std::getline(check1, intermediate, ' '))
     {
         tokens.push_back(intermediate);
     }
     int gsize= tokens.size();
     // Printing the token vector
     //for(j = 0; j < tokens.size(); j++) {
     //   std::cout << tokens[j] << '\n';
     //}
     //std::cout<< "gsize=" << gsize << "\n";
     //Graph *gc = new Graph();
     //gc->createGraph(gsize);
     int e=0;
     char str[gsize*gsize+1];
     for (j=0; j < gsize; j++) {
        for (k=0; k < gsize; k++) {
           //std::cout << "tokens[j]=" << tokens[j] << "  tokens[k]=" << tokens[k] << "\n";
           //std::cout << "maps[tokens[j]]=" << maps[tokens[j]] << "  maps[tokens[k]]=" << maps[tokens[k]] << "\n";
           if (g->hasEdge(maps[tokens[j]], maps[tokens[k]])) {
            str[j*gsize+k]='1';
            //gc->addEdge(j,k);
            e++;
           }
           else {
             str[j*gsize+k]='0';
           }
        }
    }
    str[gsize*gsize]='\0';
    //std::cout<< "e=" << e << "\n";
    //std::cout<< "str=" << str << "\n\n";
    if (e<3) continue;
    //std::cout<< "e=" << e << "\n";
    if (gsize<=16) {
        Subgraph *sg = new Subgraph(gsize);
        sg->canstr=new char[gsize*gsize+1];
        int map2[gsize];

       GraphIsomor::canonicalOrderWithMap(str, sg->canstr, gsize,map2);
       //std::cout<< "canstr=" << sg->canstr << "\n";
       //printf("\n%d\n\n",isomorSG.checkString(sg->canstr));
       if (isomorSG.checkString(sg->canstr) != 1) {
            //std::cout<< "canstr=" << sg->canstr << "\n";
            isomorSG.setString(sg->canstr, 0);
            sgv.push_back(sg);
       }
    }
    //std::cout << "sgv.size() = " << sgv.size() << std::endl;

  }
  std::vector<Subgraph *>::iterator ii;
  //std::cout << sgv.size()<< std::endl;
  for(ii=sgv.begin(); ii!=sgv.end(); ii++) {
	//std::cout << "complex=" << (*ii)->canstr << std::endl;
	Graph *gc = new Graph();
	int nswap = 0;
	int sz=sqrt(strlen((*ii)->canstr));
    Graph::strToGraph(gc, (*ii)->canstr, sz);
    int flag=0;
    //printf("no. of edges=%d\n",gc->numberEdges());
    int rn=0.25*gc->numberEdges();
    //printf("no. of edges reduction=%d\n",rn);
    int i=0, j=1;
    while(rn>0) {
        if(gc->hasEdge(i,j)) {
            gc->rmEdge(i,j);
        }
        if(j<sz-1) {
            j++;
        }
        else if(i<sz-2) {
            i++;
            j=i+1;
        }
        else flag=10;
        rn--;
    }
    do {
        Graph::randomGraph(gc, num_exchanges, num_tries, &nswap);
        char str[sz*sz+1];
        Graph::graphToStr(gc, str, sz);
        //std::cout << "random str" << i << "  " << j << " = " <<str << std::endl;
        Subgraph *sg = new Subgraph(sz);
        sg->canstr=new char[sz*sz+1];
        int map2[sz];
        GraphIsomor::canonicalOrderWithMap(str, sg->canstr, sz,map2);
        //std::cout<< "random canstr of size " << sz << " = " << sg->canstr << "\n";
        if (isomorSG.checkString(sg->canstr) != 1) {
            //std::cout<< " random generated " << "\n\n";
            isomorSG.setString(sg->canstr, 0);
            sgv_random.push_back(sg);
            flag++;
        }
        else {
            //std::cout<< " present in complex list " << "\n\n";
        }

        if(gc->hasEdge(i,j)) {
            gc->rmEdge(i,j);
        }
        if(j<sz-1) {
            j++;
        }
        else if(i<sz-2) {
            i++;
            j=i+1;
        }
        else flag=10;
        //std::cout << "sgv_random.size() = " << sgv_random.size() << std::endl;
    }while(flag<5);
    gc->deleteGraph();
    if(sgv_random.size() >= 10*sgv.size()) break;
  }

  f.close();
  mfile.close();
}

// function that returns correlation coefficient.
float correlationCoefficient(float X[], float Y[], int n)
{
    float sum_X = 0, sum_Y = 0, sum_XY = 0;
    float squareSum_X = 0, squareSum_Y = 0;

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];

        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];

        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X + X[i] * X[i];
        squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }
    if(sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y)) == 0)
        return 1;
    // use formula for calculating correlation coefficient.
    float corr = (float)(n * sum_XY - sum_X * sum_Y)
                  / sqrt((n * squareSum_X - sum_X * sum_X)
                      * (n * squareSum_Y - sum_Y * sum_Y));
    if (corr < 0)
        return -corr;
    else
        return corr;
}

void extractFeature(Graph* gc,double feature[],double norm_feature[])
{
      int train_feature[22];
      int i,j,k,l;
      int n=gc->numberNodes();
      int e=gc->numberEdges();
      //printf("n=%d   e=%d\n",n,e);


      feature[0]=n;
      feature[1]=(2.0*e)/(n*(n-1));
      norm_feature[0]=n/n;
      norm_feature[1]=(2.0*e)/(n*(n-1));
      //printf("norm_feature[0]=%lf   norm_feature[1]=%lf\n",norm_feature[0],norm_feature[1]);

      bool **adjMat = gc->adjacencyMatrix();

/**************************************************/
      double meanDegree=0.0;
      int *numnei = gc->numberNeighbours();
      for (i=0; i < n; i++) {
        meanDegree += numnei[i];
      }
      meanDegree = meanDegree/n;
      feature[2]=meanDegree;
      norm_feature[2]=meanDegree/n;
      train_feature[2]=catagorize(norm_feature[2]);
/**************************************************/

/**************************************************/
      double varDegree=0.0;
      for (i=0; i < n; i++) {
        varDegree += (numnei[i]-meanDegree)*(numnei[i]-meanDegree);
      }
      varDegree=varDegree/n;
      feature[3]=varDegree;
      norm_feature[3]=varDegree/n;
      train_feature[3]=catagorize(norm_feature[3]);
/**************************************************/

/**************************************************/
      int degrees[n];
      for (i=0; i < n; i++) {
        degrees[i] = numnei[i];
      }
      for (i=0; i < n-1; i++) {
            for (j=n-1; j > i; j--) {
                if(degrees[j] > degrees[j-1]) {
                    int temp=degrees[j];
                    degrees[j]=degrees[j-1];
                    degrees[j-1]=temp;
                }
            }
      }
      /*for (i=0; i < n; i++) {
        std::cout << degrees[i] << "  ";
      }
      std::cout << std::endl;*/
      int medianDegree=degrees[n/2];
      feature[4]=medianDegree;
      norm_feature[4]=1.0*medianDegree/n;
      train_feature[4]=catagorize(norm_feature[4]);
/**************************************************/

/**************************************************/
      int maxDegree=numnei[0];
      for (i=1; i < n; i++) {
        if (maxDegree < numnei[i]) {
            maxDegree = numnei[i];
        }
      }
      feature[5]=maxDegree;
      norm_feature[5]=1.0*maxDegree/n;
      train_feature[5]=catagorize(norm_feature[5]);
/**************************************************/

/**************************************************/
      int xx[2*e], yy[2*e], kk=0;
      double corrxy[2*e], corrx[2*e], corry[2*e];
      double meanxx=0.0, meanyy=0.0;
      for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                    if(adjMat[i][j]) {
                        xx[kk]=numnei[i];
                        meanxx+=xx[kk];
                        yy[kk]=numnei[j];
                        kk++;
                    }
            }
      }
      //if (kk!=2*e) {printf("\nerror\n");}
      //else printf("kk=%d\n",kk);
      meanxx= meanxx/kk;
      meanyy=meanxx;
      /*std::cout << std::endl << std::endl;
      for(i=0;i<kk;i++) {
            std::cout << xx[i] << "   " << yy[i] << "\n";
      }
      std::cout << std::endl << std::endl;
      std::cout << "meanx=" << meanxx << "    meany=" << meanyy << "\n";*/

      double maxDegreeCorrelation=-5.0;
      double meanDegreeCorrelation=0.0;
      double corrxsum=0.0;
      double corrysum=0.0;
      double varDegreeCorrelation;
      for (i=0; i < kk; i++) {
        corrxy[i]=(xx[i]-meanxx)*(yy[i]-meanyy);
        corrx[i]=(xx[i]-meanxx)*(xx[i]-meanxx);
        corry[i]=(yy[i]-meanyy)*(yy[i]-meanyy);
        if(corrxy[i] > maxDegreeCorrelation) maxDegreeCorrelation=corrxy[i];
        meanDegreeCorrelation+=corrxy[i];
        corrxsum+=corrx[i];
        corrysum+=corry[i];
      }
      if (meanDegreeCorrelation==0) varDegreeCorrelation=0;
      else varDegreeCorrelation=meanDegreeCorrelation/sqrt(corrxsum*corrysum);
      meanDegreeCorrelation=meanDegreeCorrelation/kk;


      /*std::cout << std::endl << std::endl;
      for(i=0;i<kk;i++) {
            std::cout << "corrx=" << corrx[i] << "   " << "corry=" << corry[i] << "   " << "corrxy=" << corrxy[i] << "\n";
      }
      std::cout << std::endl << std::endl;*/
      feature[6]=meanDegreeCorrelation;
      feature[7]=varDegreeCorrelation;
      feature[8]=maxDegreeCorrelation;
      norm_feature[6]=meanDegreeCorrelation/n;
      norm_feature[7]=varDegreeCorrelation/n;
      norm_feature[8]=maxDegreeCorrelation/n;
/**************************************************/

      double cc[n];
      for (i=0; i < n; i++) {
        if(numnei[i]==0 || numnei[i]==1) cc[i]=0;
        else {
        int ne=0;
        int *temp = gc->getNeighboursArray(i);
        for (j=0; j < numnei[i]-1; j++) {
                for (k=j; k < numnei[i]; k++) {
                        if(gc->hasEdge(temp[j],temp[k])) ne++;
                }
        }
        cc[i] = (2.0*ne)/(numnei[i]*(numnei[i]-1));
        }
      }
      /*for(i=0;i<n;i++) {
            std::cout << cc[i] << "  ";
      }
      std::cout << std::endl;*/

      double meanClusteringCoeff=0.0;
      for (i=0; i < n; i++) {
        meanClusteringCoeff += cc[i];
      }
      meanClusteringCoeff = meanClusteringCoeff/n;

      double varClusteringCoeff=0.0;
      for (i=0; i < n; i++) {
        varClusteringCoeff += (cc[i]-meanClusteringCoeff)*(cc[i]-meanClusteringCoeff);
      }
      varClusteringCoeff=varClusteringCoeff/n;

      double maxClusteringCoeff=cc[0];
      for (i=1; i < n; i++) {
        if (maxClusteringCoeff < cc[i])
            maxClusteringCoeff = cc[i];
      }

      feature[9]=meanClusteringCoeff;
      feature[10]=varClusteringCoeff;
      feature[11]=maxClusteringCoeff;
      norm_feature[9]=meanClusteringCoeff;
      norm_feature[10]=varClusteringCoeff;
      norm_feature[11]=maxClusteringCoeff;
/**************************************************/


      int cn[n][n]; //common neighbour
      for (i=0; i < n; i++) {
            for (j=0; j < n; j++) {
                cn[i][j]=0;
                if(i==j) continue;
                for (k=0; k < n; k++) {
                    if (k==i || k==j) continue;
                    if(gc->hasEdge(i,k) && gc->hasEdge(j,k)) cn[i][j]++;
                }
            }
      }
      /*for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                    std::cout << cn[i][j] << "  ";
            }
            std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;*/
      double avgcn[n]; //average common neighbour
      for (i=0; i < n; i++) {
            avgcn[i]=0.0;
            for (j=0; j < n; j++) {
                avgcn[i]+=cn[i][j];
            }
            avgcn[i]=avgcn[i]/n;
      }
      double tc[n];
      for (i=0; i < n; i++) {
        if(numnei[i]==0 || numnei[i]==1) tc[i]=0;
        else tc[i]=avgcn[i]/numnei[i];
      }
      /*for(i=0;i<n;i++) {
            std::cout << tc[i] << "  ";
      }
      std::cout << std::endl;*/
      double meanTopologicCoeff=0.0;
      for (i=0; i < n; i++) {
        meanTopologicCoeff += tc[i];
      }
      meanTopologicCoeff = meanTopologicCoeff/n;

      double varTopologicCoeff=0.0;
      for (i=0; i < n; i++) {
        varTopologicCoeff += (tc[i]-meanTopologicCoeff)*(tc[i]-meanTopologicCoeff);
      }
      varTopologicCoeff=varTopologicCoeff/n;

      double maxTopologicCoeff=tc[0];
      for (i=1; i < n; i++) {
        if (maxTopologicCoeff < tc[i])
            maxTopologicCoeff = tc[i];
      }

      feature[12]=meanTopologicCoeff;
      feature[13]=varTopologicCoeff;
      feature[14]=maxTopologicCoeff;
      norm_feature[12]=meanTopologicCoeff;
      norm_feature[13]=varTopologicCoeff;
      norm_feature[14]=maxTopologicCoeff;
/**************************************************/


      double A[n][n];
      double U[n][n];
      double V[n][n];
      double singular_values[n];
      double* dummy_array;

      //(your code to initialize the matrix A)
      for (i=0; i < n; i++) {
            for (j=0; j < n; j++) {
                A[i][j]=adjMat[i][j];
            }
      }
      dummy_array = (double*) malloc(n * sizeof(double));
      if (dummy_array == NULL) {printf(" No memory available\n"); exit(0); }

      double err = Singular_Value_Decomposition((double*) A,n,n,(double*)U,singular_values,(double*)V,dummy_array);

      free(dummy_array);
      if (err < 0) {
            //printf(" Failed to converge\n");
            feature[15]=0;
            feature[16]=0;
            feature[17]=0;
            norm_feature[15]=0;
            norm_feature[16]=0;
            norm_feature[17]=0;
      }
      else {
            /*printf(" The singular value decomposition of A is \n");
            for(i=0;i<n;i++) {
                std::cout << singular_values[i] << "  ";
            }
            std::cout << std::endl;*/
            double eigenValue_1=singular_values[0];
            double eigenValue_2=singular_values[1];
            double eigenValue_3=singular_values[2];
            feature[15]=eigenValue_1;
            feature[16]=eigenValue_2;
            feature[17]=eigenValue_3;
            norm_feature[15]=eigenValue_1/n;
            norm_feature[16]=eigenValue_2/n;
            norm_feature[17]=eigenValue_3/n;
      }
  /**************************************************/


      int matrix[30][30];
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (adjMat[i][j] == 0) matrix[i][j]=999;
            else if (i==j) matrix[i][j]=0;
            else matrix[i][j]=adjMat[i][j];
        }
      }
      floydWarshall(matrix, n);
      double maxLength=0.0;
      int dc=0;
      double aveLength=0.0;
      for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            if (matrix[i][j]==999) continue;
            aveLength=aveLength+matrix[i][j];
            dc++;
            //printf("dc=%d   d=%d   aveL=%lf\n", dc,matrix[i][j],aveLength);
            if (matrix[i][j] > maxLength) maxLength=matrix[i][j];
        }
      }
      aveLength=aveLength/dc;
      //printf("%maxLength=%d\n", maxLength);
      //printf("aveLength=%lf\n", aveLength);
      feature[18]=maxLength;
      feature[19]=aveLength;
      norm_feature[18]=maxLength/(n-1);
      norm_feature[19]=aveLength/(n-1);
  /**************************************************/


/*For an input graph G = (V, E), we assign the weight of an edge [u, v]
to be the number of neighbors shared by the vertices u and v.
We define the weight of each vertex
to be the sum of the weights of its incident edges.*/
      int W[n];
      double maxWeight=0.0;
      double aveWeight=0.0;
      for(i=0;i<n;i++) {
        W[i]=0;
        int *temp = gc->getNeighboursArray(i);
        //for(k=0;k<numnei[i];k++) printf("%d  ",temp[k]);
        //printf("\n\n");
        for(j=0;j<numnei[i];j++) {
            int *temp2 = gc->getNeighboursArray(temp[j]);
            //for(k=0;k<numnei[temp[j]];k++) printf("%d  ",temp2[k]);
            //printf("\n");
            for(int k1=0;k1<numnei[i];k1++) {
                for(int k2=0;k2<numnei[temp[j]];k2++) {
                    if(temp[k1]==temp2[k2]) W[i]++;
                }
            }
        }
        //printf("W=%d\n",W[i]);
        //printf("\n\n");
      }
      for(i=0;i<n;i++) {
            if(W[i]>maxWeight) maxWeight=W[i];
            aveWeight=aveWeight+W[i];
            //printf("%d  ",W[i]);
      }
      //printf("\n");
      aveWeight=aveWeight/n;
      feature[20]=maxWeight;
      feature[21]=aveWeight;
      norm_feature[20]=maxWeight/((n-1)*(n-2));
      norm_feature[21]=aveWeight/((n-1)*(n-2));
//      for(i=0;i<22;i++) {
//        std::cout<< std::setw(12) << norm_feature[i] << "  ";
//      }
}

// Implementing floyd warshall algorithm
void floydWarshall(int matrix[][30], int nV) {
  int i, j, k;
  for (k = 0; k < nV; k++) {
    for (i = 0; i < nV; i++) {
      for (j = 0; j < nV; j++) {
        if (matrix[i][k] + matrix[k][j] < matrix[i][j])
          matrix[i][j] = matrix[i][k] + matrix[k][j];
      }
    }
  }

}
int catagorize(double v)
{
    if(v<= 0.01) return 0;
    else if(v <=0.1) return 1;
    else if(v <=0.2) return 2;
    else if(v <=0.3) return 3;
    else if(v <=0.4) return 4;
    else if(v <=0.5) return 5;
    else if(v <=0.6) return 6;
    else if(v <=0.7) return 7;
    else if(v <=0.8) return 8;
    else if(v <=0.9) return 9;
    else if(v <=1) return 10;
}

////////////////////////////////////////////////////////////////////////////////
//  int Singular_Value_Decomposition(double* A, int nrows, int ncols,         //
//        double* U, double* singular_values, double* V, double* dummy_array) //
//                                                                            //
//  Description:                                                              //
//     This routine decomposes an m x n matrix A, with m >= n, into a product //
//     of the three matrices U, D, and V', i.e. A = UDV', where U is an m x n //
//     matrix whose columns are orthogonal, D is a n x n diagonal matrix, and //
//     V is an n x n orthogonal matrix.  V' denotes the transpose of V.  If   //
//     m < n, then the procedure may be used for the matrix A'.  The singular //
//     values of A are the diagonal elements of the diagonal matrix D and     //
//     correspond to the positive square roots of the eigenvalues of the      //
//     matrix A'A.                                                            //
//                                                                            //
//     This procedure programmed here is based on the method of Golub and     //
//     Reinsch as given on pages 134 - 151 of the "Handbook for Automatic     //
//     Computation vol II - Linear Algebra" edited by Wilkinson and Reinsch   //
//     and published by Springer-Verlag, 1971.                                //
//                                                                            //
//     The Golub and Reinsch's method for decomposing the matrix A into the   //
//     product U, D, and V' is performed in three stages:                     //
//       Stage 1:  Decompose A into the product of three matrices U1, B, V1'  //
//         A = U1 B V1' where B is a bidiagonal matrix, and U1, and V1 are a  //
//         product of Householder transformations.                            //
//       Stage 2:  Use Given' transformations to reduce the bidiagonal matrix //
//         B into the product of the three matrices U2, D, V2'.  The singular //
//         value decomposition is then UDV'where U = U2 U1 and V' = V1' V2'.  //
//       Stage 3:  Sort the matrix D in decreasing order of the singular      //
//         values and interchange the columns of both U and V to reflect any  //
//         change in the order of the singular values.                        //
//                                                                            //
//     After performing the singular value decomposition for A, call          //
//     Singular_Value_Decomposition to solve the equation Ax = B or call      //
//     Singular_Value_Decomposition_Inverse to calculate the pseudo-inverse   //
//     of A.                                                                  //
//                                                                            //
//  Arguments:                                                                //
//     double* A                                                              //
//        On input, the pointer to the first element of the matrix            //
//        A[nrows][ncols].  The matrix A is unchanged.                        //
//     int nrows                                                              //
//        The number of rows of the matrix A.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix A.                              //
//     double* U                                                              //
//        On input, a pointer to a matrix with the same number of rows and    //
//        columns as the matrix A.  On output, the matrix with mutually       //
//        orthogonal columns which is the left-most factor in the singular    //
//        value decomposition of A.                                           //
//     double* singular_values                                                //
//        On input, a pointer to an array dimensioned to same as the number   //
//        of columns of the matrix A, ncols.  On output, the singular values  //
//        of the matrix A sorted in decreasing order.  This array corresponds //
//        to the diagonal matrix in the singular value decomposition of A.    //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix A, i.e. V[ncols][ncols].   //
//        On output, the orthogonal matrix whose transpose is the right-most  //
//        factor in the singular value decomposition of A.                    //
//     double* dummy_array                                                    //
//        On input, a pointer to an array dimensioned to same as the number   //
//        of columns of the matrix A, ncols.  This array is used to store     //
//        the super-diagonal elements resulting from the Householder reduction//
//        of the matrix A to bidiagonal form.  And as an input to the Given's //
//        procedure to reduce the bidiagonal form to diagonal form.           //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - During the Given's reduction of the bidiagonal form to    //
//                  diagonal form the procedure failed to terminate within    //
//                  MAX_ITERATION_COUNT iterations.                           //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double A[M][N];                                                        //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double singular_values[N];                                             //
//     double* dummy_array;                                                   //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//     dummy_array = (double*) malloc(N * sizeof(double));                    //
//     if (dummy_array == NULL) {printf(" No memory available\n"); exit(0); } //
//                                                                            //
//     err = Singular_Value_Decomposition((double*) A, M, N, (double*) U,     //
//                              singular_values, (double*) V, dummy_array);   //
//                                                                            //
//     free(dummy_array);                                                     //
//     if (err < 0) printf(" Failed to converge\n");                          //
//     else { printf(" The singular value decomposition of A is \n");         //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Singular_Value_Decomposition(double* A, int nrows, int ncols, double* U,
                      double* singular_values, double* V, double* dummy_array)
{
   Householders_Reduction_to_Bidiagonal_Form( A, nrows, ncols, U, V,
                                                singular_values, dummy_array);

   if (Givens_Reduction_to_Diagonal_Form( nrows, ncols, U, V,
                                singular_values, dummy_array ) < 0) return -1;

   Sort_by_Decreasing_Singular_Values(nrows, ncols, singular_values, U, V);

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
// static void Householders_Reduction_to_Bidiagonal_Form(double* A, int nrows,//
//  int ncols, double* U, double* V, double* diagonal, double* superdiagonal )//
//                                                                            //
//  Description:                                                              //
//     This routine decomposes an m x n matrix A, with m >= n, into a product //
//     of the three matrices U, B, and V', i.e. A = UBV', where U is an m x n //
//     matrix whose columns are orthogonal, B is a n x n bidiagonal matrix,   //
//     and V is an n x n orthogonal matrix.  V' denotes the transpose of V.   //
//     If m < n, then the procedure may be used for the matrix A'.  The       //
//                                                                            //
//     The matrix U is the product of Householder transformations which       //
//     annihilate the subdiagonal components of A while the matrix V is       //
//     the product of Householder transformations which annihilate the        //
//     components of A to the right of the superdiagonal.                     //
//                                                                            //
//     The Householder transformation which leaves invariant the first k-1    //
//     elements of the k-th column and annihilates the all the elements below //
//     the diagonal element is P = I - (2/u'u)uu', u is an nrows-dimensional  //
//     vector the first k-1 components of which are zero and the last         //
//     components agree with the current transformed matrix below the diagonal//
//     diagonal, the remaining k-th element is the diagonal element - s, where//
//     s = (+/-)sqrt(sum of squares of the elements below the diagonal), the  //
//     sign is chosen opposite that of the diagonal element.                  //
//                                                                            //
//  Arguments:                                                                //
//     double* A                                                              //
//        On input, the pointer to the first element of the matrix            //
//        A[nrows][ncols].  The matrix A is unchanged.                        //
//     int nrows                                                              //
//        The number of rows of the matrix A.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix A.                              //
//     double* U                                                              //
//        On input, a pointer to a matrix with the same number of rows and    //
//        columns as the matrix A.  On output, the matrix with mutually       //
//        orthogonal columns which is the left-most factor in the bidiagonal  //
//        decomposition of A.                                                 //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix A, i.e. V[ncols][ncols].   //
//        On output, the orthogonal matrix whose transpose is the right-most  //
//        factor in the bidiagonal decomposition of A.                        //
//     double* diagonal                                                       //
//        On input, a pointer to an array dimensioned to same as the number   //
//        of columns of the matrix A, ncols.  On output, the diagonal of the  //
//        bidiagonal matrix.                                                  //
//     double* superdiagonal                                                  //
//        On input, a pointer to an array dimensioned to same as the number   //
//        of columns of the matrix A, ncols.  On output, the superdiagonal    //
//        of the bidiagonal matrix.                                           //
//                                                                            //
//  Return Values:                                                            //
//     The function is of type void and therefore does not return a value.    //
//     The matrices U, V, and the diagonal and superdiagonal are calculated   //
//     using the addresses passed in the argument list.                       //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double A[M][N];                                                        //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double diagonal[N];                                                    //
//     double superdiagonal[N];                                               //
//                                                                            //
//     (your code to initialize the matrix A - Note this routine is not       //
//     (accessible from outside i.e. it is declared static)                   //
//                                                                            //
//     Householders_Reduction_to_Bidiagonal_Form((double*) A, nrows, ncols,   //
//                   (double*) U, (double*) V, diagonal, superdiagonal )      //
//                                                                            //
//     free(dummy_array);                                                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Householders_Reduction_to_Bidiagonal_Form(double* A, int nrows,
    int ncols, double* U, double* V, double* diagonal, double* superdiagonal )
{
   int i,j,k,ip1;
   double s, s2, si, scale;
   double dum;
   double *pu, *pui, *pv, *pvi;
   double half_norm_squared;

// Copy A to U

   memcpy(U,A, sizeof(double) * nrows * ncols);

//

   diagonal[0] = 0.0;
   s = 0.0;
   scale = 0.0;
   for ( i = 0, pui = U, ip1 = 1; i < ncols; pui += ncols, i++, ip1++ ) {
      superdiagonal[i] = scale * s;
//
//                  Perform Householder transform on columns.
//
//       Calculate the normed squared of the i-th column vector starting at
//       row i.
//
      for (j = i, pu = pui, scale = 0.0; j < nrows; j++, pu += ncols)
         scale += fabs( *(pu + i) );

      if (scale > 0.0) {
         for (j = i, pu = pui, s2 = 0.0; j < nrows; j++, pu += ncols) {
            *(pu + i) /= scale;
            s2 += *(pu + i) * *(pu + i);
         }
//
//
//       Chose sign of s which maximizes the norm
//
         s = ( *(pui + i) < 0.0 ) ? sqrt(s2) : -sqrt(s2);
//
//       Calculate -2/u'u
//
         half_norm_squared = *(pui + i) * s - s2;
//
//       Transform remaining columns by the Householder transform.
//
         *(pui + i) -= s;

         for (j = ip1; j < ncols; j++) {
            for (k = i, si = 0.0, pu = pui; k < nrows; k++, pu += ncols)
               si += *(pu + i) * *(pu + j);
            si /= half_norm_squared;
            for (k = i, pu = pui; k < nrows; k++, pu += ncols) {
               *(pu + j) += si * *(pu + i);
            }
         }
      }
      for (j = i, pu = pui; j < nrows; j++, pu += ncols) *(pu + i) *= scale;
      diagonal[i] = s * scale;
//
//                  Perform Householder transform on rows.
//
//       Calculate the normed squared of the i-th row vector starting at
//       column i.
//
      s = 0.0;
      scale = 0.0;
      if (i >= nrows || i == (ncols - 1) ) continue;
      for (j = ip1; j < ncols; j++) scale += fabs ( *(pui + j) );
      if ( scale > 0.0 ) {
         for (j = ip1, s2 = 0.0; j < ncols; j++) {
            *(pui + j) /= scale;
            s2 += *(pui + j) * *(pui + j);
         }
         s = ( *(pui + ip1) < 0.0 ) ? sqrt(s2) : -sqrt(s2);
//
//       Calculate -2/u'u
//
         half_norm_squared = *(pui + ip1) * s - s2;
//
//       Transform the rows by the Householder transform.
//
         *(pui + ip1) -= s;
         for (k = ip1; k < ncols; k++)
            superdiagonal[k] = *(pui + k) / half_norm_squared;
         if ( i < (nrows - 1) ) {
            for (j = ip1, pu = pui + ncols; j < nrows; j++, pu += ncols) {
               for (k = ip1, si = 0.0; k < ncols; k++)
                  si += *(pui + k) * *(pu + k);
               for (k = ip1; k < ncols; k++) {
                  *(pu + k) += si * superdiagonal[k];
               }
            }
         }
         for (k = ip1; k < ncols; k++) *(pui + k) *= scale;
      }
   }

// Update V
   pui = U + ncols * (ncols - 2);
   pvi = V + ncols * (ncols - 1);
   *(pvi + ncols - 1) = 1.0;
   s = superdiagonal[ncols - 1];
   pvi -= ncols;
   for (i = ncols - 2, ip1 = ncols - 1; i >= 0; i--, pui -= ncols,
                                                      pvi -= ncols, ip1-- ) {
      if ( s != 0.0 ) {
         pv = pvi + ncols;
         for (j = ip1; j < ncols; j++, pv += ncols)
            *(pv + i) = ( *(pui + j) / *(pui + ip1) ) / s;
         for (j = ip1; j < ncols; j++) {
            si = 0.0;
            for (k = ip1, pv = pvi + ncols; k < ncols; k++, pv += ncols)
               si += *(pui + k) * *(pv + j);
            for (k = ip1, pv = pvi + ncols; k < ncols; k++, pv += ncols)
               *(pv + j) += si * *(pv + i);
         }
      }
      pv = pvi + ncols;
      for ( j = ip1; j < ncols; j++, pv += ncols ) {
         *(pvi + j) = 0.0;
         *(pv + i) = 0.0;
      }
      *(pvi + i) = 1.0;
      s = superdiagonal[i];
   }

// Update U

   pui = U + ncols * (ncols - 1);
   for (i = ncols - 1, ip1 = ncols; i >= 0; ip1 = i, i--, pui -= ncols ) {
      s = diagonal[i];
      for ( j = ip1; j < ncols; j++) *(pui + j) = 0.0;
      if ( s != 0.0 ) {
         for (j = ip1; j < ncols; j++) {
            si = 0.0;
            pu = pui + ncols;
            for (k = ip1; k < nrows; k++, pu += ncols)
               si += *(pu + i) * *(pu + j);
            si = (si / *(pui + i) ) / s;
            for (k = i, pu = pui; k < nrows; k++, pu += ncols)
               *(pu + j) += si * *(pu + i);
         }
         for (j = i, pu = pui; j < nrows; j++, pu += ncols){
            *(pu + i) /= s;
         }
      }
      else
         for (j = i, pu = pui; j < nrows; j++, pu += ncols) *(pu + i) = 0.0;
      *(pui + i) += 1.0;
   }
}


////////////////////////////////////////////////////////////////////////////////
// static int Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,        //
//         double* U, double* V, double* diagonal, double* superdiagonal )    //
//                                                                            //
//  Description:                                                              //
//     This routine decomposes a bidiagonal matrix given by the arrays        //
//     diagonal and superdiagonal into a product of three matrices U1, D and  //
//     V1', the matrix U1 premultiplies U and is returned in U, the matrix    //
//     V1 premultiplies V and is returned in V.  The matrix D is a diagonal   //
//     matrix and replaces the array diagonal.                                //
//                                                                            //
//     The method used to annihilate the offdiagonal elements is a variant    //
//     of the QR transformation.  The method consists of applying Givens      //
//     rotations to the right and the left of the current matrix until        //
//     the new off-diagonal elements are chased out of the matrix.            //
//                                                                            //
//     The process is an iterative process which due to roundoff errors may   //
//     not converge within a predefined number of iterations.  (This should   //
//     be unusual.)                                                           //
//                                                                            //
//  Arguments:                                                                //
//     int nrows                                                              //
//        The number of rows of the matrix U.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix U.                              //
//     double* U                                                              //
//        On input, a pointer to a matrix already initialized to a matrix     //
//        with mutually orthogonal columns.   On output, the matrix with      //
//        mutually orthogonal columns.                                        //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix U, i.e. V[ncols][ncols].   //
//        The matrix V is assumed to be initialized to an orthogonal matrix.  //
//        On output, V is an orthogonal matrix.                               //
//     double* diagonal                                                       //
//        On input, a pointer to an array of dimension ncols which initially  //
//        contains the diagonal of the bidiagonal matrix.  On output, the     //
//        it contains the diagonal of the diagonal matrix.                    //
//     double* superdiagonal                                                  //
//        On input, a pointer to an array of dimension ncols which initially  //
//        the first component is zero and the successive components form the  //
//        superdiagonal of the bidiagonal matrix.                             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The procedure failed to terminate within                  //
//                  MAX_ITERATION_COUNT iterations.                           //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double diagonal[N];                                                    //
//     double superdiagonal[N];                                               //
//     int err;                                                               //
//                                                                            //
//     (your code to initialize the matrices U, V, diagonal, and )            //
//     ( superdiagonal.  - Note this routine is not accessible from outside)  //
//     ( i.e. it is declared static.)                                         //
//                                                                            //
//     err = Givens_Reduction_to_Diagonal_Form( M,N,(double*)U,(double*)V,    //
//                                                 diagonal, superdiagonal ); //
//     if ( err < 0 ) printf("Failed to converge\n");                         //
//     else { ... }                                                           //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static int Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,
           double* U, double* V, double* diagonal, double* superdiagonal )
{

   double epsilon;
   double c, s;
   double f,g,h;
   double x,y,z;
   double *pu, *pv;
   int i,j,k,m;
   int rotation_test;
   int iteration_count;

   for (i = 0, x = 0.0; i < ncols; i++) {
      y = fabs(diagonal[i]) + fabs(superdiagonal[i]);
      if ( x < y ) x = y;
   }
   epsilon = x * DBL_EPSILON;
   for (k = ncols - 1; k >= 0; k--) {
      iteration_count = 0;
      while(1) {
         rotation_test = 1;
         for (m = k; m >= 0; m--) {
            if (fabs(superdiagonal[m]) <= epsilon) {rotation_test = 0; break;}
            if (fabs(diagonal[m-1]) <= epsilon) break;
         }
         if (rotation_test) {
            c = 0.0;
            s = 1.0;
            for (i = m; i <= k; i++) {
               f = s * superdiagonal[i];
               superdiagonal[i] *= c;
               if (fabs(f) <= epsilon) break;
               g = diagonal[i];
               h = sqrt(f*f + g*g);
               diagonal[i] = h;
               c = g / h;
               s = -f / h;
               for (j = 0, pu = U; j < nrows; j++, pu += ncols) {
                  y = *(pu + m - 1);
                  z = *(pu + i);
                  *(pu + m - 1 ) = y * c + z * s;
                  *(pu + i) = -y * s + z * c;
               }
            }
         }
         z = diagonal[k];
         if (m == k ) {
            if ( z < 0.0 ) {
               diagonal[k] = -z;
               for ( j = 0, pv = V; j < ncols; j++, pv += ncols)
                  *(pv + k) = - *(pv + k);
            }
            break;
         }
         else {
            if ( iteration_count >= MAX_ITERATION_COUNT ) return -1;
            iteration_count++;
            x = diagonal[m];
            y = diagonal[k-1];
            g = superdiagonal[k-1];
            h = superdiagonal[k];
            f = ( (y - z) * ( y + z ) + (g - h) * (g + h) )/(2.0 * h * y);
            g = sqrt( f * f + 1.0 );
            if ( f < 0.0 ) g = -g;
            f = ( (x - z) * (x + z) + h * (y / (f + g) - h) ) / x;
// Next QR Transformtion
            c = 1.0;
            s = 1.0;
            for (i = m + 1; i <= k; i++) {
               g = superdiagonal[i];
               y = diagonal[i];
               h = s * g;
               g *= c;
               z = sqrt( f * f + h * h );
               superdiagonal[i-1] = z;
               c = f / z;
               s = h / z;
               f =  x * c + g * s;
               g = -x * s + g * c;
               h = y * s;
               y *= c;
               for (j = 0, pv = V; j < ncols; j++, pv += ncols) {
                  x = *(pv + i - 1);
                  z = *(pv + i);
                  *(pv + i - 1) = x * c + z * s;
                  *(pv + i) = -x * s + z * c;
               }
               z = sqrt( f * f + h * h );
               diagonal[i - 1] = z;
               if (z != 0.0) {
                  c = f / z;
                  s = h / z;
               }
               f = c * g + s * y;
               x = -s * g + c * y;
               for (j = 0, pu = U; j < nrows; j++, pu += ncols) {
                  y = *(pu + i - 1);
                  z = *(pu + i);
                  *(pu + i - 1) = c * y + s * z;
                  *(pu + i) = -s * y + c * z;
               }
            }
            superdiagonal[m] = 0.0;
            superdiagonal[k] = f;
            diagonal[k] = x;
         }
      }
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
// static void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,       //
//                            double* singular_values, double* U, double* V)  //
//                                                                            //
//  Description:                                                              //
//     This routine sorts the singular values from largest to smallest        //
//     singular value and interchanges the columns of U and the columns of V  //
//     whenever a swap is made.  I.e. if the i-th singular value is swapped   //
//     with the j-th singular value, then the i-th and j-th columns of U are  //
//     interchanged and the i-th and j-th columns of V are interchanged.      //
//                                                                            //
//  Arguments:                                                                //
//     int nrows                                                              //
//        The number of rows of the matrix U.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix U.                              //
//     double* singular_values                                                //
//        On input, a pointer to the array of singular values.  On output, the//
//        sorted array of singular values.                                    //
//     double* U                                                              //
//        On input, a pointer to a matrix already initialized to a matrix     //
//        with mutually orthogonal columns.  On output, the matrix with       //
//        mutually orthogonal possibly permuted columns.                      //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix U, i.e. V[ncols][ncols].   //
//        The matrix V is assumed to be initialized to an orthogonal matrix.  //
//        On output, V is an orthogonal matrix with possibly permuted columns.//
//                                                                            //
//  Return Values:                                                            //
//        The function is of type void.                                       //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double diagonal[N];                                                    //
//                                                                            //
//     (your code to initialize the matrices U, V, and diagonal. )            //
//     ( - Note this routine is not accessible from outside)                  //
//     ( i.e. it is declared static.)                                         //
//                                                                            //
//     Sort_by_Decreasing_Singular_Values(nrows, ncols, singular_values,      //
//                                                 (double*) U, (double*) V); //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,
                                double* singular_values, double* U, double* V)
{
   int i,j,max_index;
   double temp;
   double *p1, *p2;

   for (i = 0; i < ncols - 1; i++) {
      max_index = i;
      for (j = i + 1; j < ncols; j++)
         if (singular_values[j] > singular_values[max_index] )
            max_index = j;
      if (max_index == i) continue;
      temp = singular_values[i];
      singular_values[i] = singular_values[max_index];
      singular_values[max_index] = temp;
      p1 = U + max_index;
      p2 = U + i;
      for (j = 0; j < nrows; j++, p1 += ncols, p2 += ncols) {
         temp = *p1;
         *p1 = *p2;
         *p2 = temp;
      }
      p1 = V + max_index;
      p2 = V + i;
      for (j = 0; j < ncols; j++, p1 += ncols, p2 += ncols) {
         temp = *p1;
         *p1 = *p2;
         *p2 = temp;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
//  void Singular_Value_Decomposition_Solve(double* U, double* D, double* V,  //
//              double tolerance, int nrows, int ncols, double *B, double* x) //
//                                                                            //
//  Description:                                                              //
//     This routine solves the system of linear equations Ax=B where A =UDV', //
//     is the singular value decomposition of A.  Given UDV'x=B, then         //
//     x = V(1/D)U'B, where 1/D is the pseudo-inverse of D, i.e. if D[i] > 0  //
//     then (1/D)[i] = 1/D[i] and if D[i] = 0, then (1/D)[i] = 0.  Since      //
//     the singular values are subject to round-off error.  A tolerance is    //
//     given so that if D[i] < tolerance, D[i] is treated as if it is 0.      //
//     The default tolerance is D[0] * DBL_EPSILON * ncols, if the user       //
//     specified tolerance is less than the default tolerance, the default    //
//     tolerance is used.                                                     //
//                                                                            //
//  Arguments:                                                                //
//     double* U                                                              //
//        A matrix with mutually orthonormal columns.                         //
//     double* D                                                              //
//        A diagonal matrix with decreasing non-negative diagonal elements.   //
//        i.e. D[i] > D[j] if i < j and D[i] >= 0 for all i.                  //
//     double* V                                                              //
//        An orthogonal matrix.                                               //
//     double tolerance                                                       //
//        An lower bound for non-zero singular values (provided tolerance >   //
//        ncols * DBL_EPSILON * D[0]).                                        //
//     int nrows                                                              //
//        The number of rows of the matrix U and B.                           //
//     int ncols                                                              //
//        The number of columns of the matrix U.  Also the number of rows and //
//        columns of the matrices D and V.                                    //
//     double* B                                                              //
//        A pointer to a vector dimensioned as nrows which is the  right-hand //
//        side of the equation Ax = B where A = UDV'.                         //
//     double* x                                                              //
//        A pointer to a vector dimensioned as ncols, which is the least      //
//        squares solution of the equation Ax = B where A = UDV'.             //
//                                                                            //
//  Return Values:                                                            //
//        The function is of type void.                                       //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     #define NB                                                             //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double D[N];                                                           //
//     double B[M];                                                           //
//     double x[N];                                                           //
//     double tolerance;                                                      //
//                                                                            //
//     (your code to initialize the matrices U,D,V,B)                         //
//                                                                            //
//     Singular_Value_Decomposition_Solve((double*) U, D, (double*) V,        //
//                                              tolerance, M, N, B, x, bcols) //
//                                                                            //
//     printf(" The solution of Ax=B is \n");                                 //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //

void Singular_Value_Decomposition_Solve(double* U, double* D, double* V,
                double tolerance, int nrows, int ncols, double *B, double* x)
{
   int i,j,k;
   double *pu, *pv;
   double dum;

   dum = DBL_EPSILON * D[0] * (double) ncols;
   if (tolerance < dum) tolerance = dum;

   for ( i = 0, pv = V; i < ncols; i++, pv += ncols) {
      x[i] = 0.0;
      for (j = 0; j < ncols; j++)
         if (D[j] > tolerance ) {
            for (k = 0, dum = 0.0, pu = U; k < nrows; k++, pu += ncols)
               dum += *(pu + j) * B[k];
            x[i] += dum * *(pv + j) / D[j];
         }
   }
}


////////////////////////////////////////////////////////////////////////////////
//  void Singular_Value_Decomposition_Inverse(double* U, double* D, double* V,//
//                     double tolerance, int nrows, int ncols, double *Astar) //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the pseudo-inverse of the matrix A = UDV'.     //
//     where U, D, V constitute the singular value decomposition of A.        //
//     Let Astar be the pseudo-inverse then Astar = V(1/D)U', where 1/D is    //
//     the pseudo-inverse of D, i.e. if D[i] > 0 then (1/D)[i] = 1/D[i] and   //
//     if D[i] = 0, then (1/D)[i] = 0.  Because the singular values are       //
//     subject to round-off error.  A tolerance is given so that if           //
//     D[i] < tolerance, D[i] is treated as if it were 0.                     //
//     The default tolerance is D[0] * DBL_EPSILON * ncols, assuming that the //
//     diagonal matrix of singular values is sorted from largest to smallest, //
//     if the user specified tolerance is less than the default tolerance,    //
//     then the default tolerance is used.                                    //
//                                                                            //
//  Arguments:                                                                //
//     double* U                                                              //
//        A matrix with mutually orthonormal columns.                         //
//     double* D                                                              //
//        A diagonal matrix with decreasing non-negative diagonal elements.   //
//        i.e. D[i] > D[j] if i < j and D[i] >= 0 for all i.                  //
//     double* V                                                              //
//        An orthogonal matrix.                                               //
//     double tolerance                                                       //
//        An lower bound for non-zero singular values (provided tolerance >   //
//        ncols * DBL_EPSILON * D[0]).                                        //
//     int nrows                                                              //
//        The number of rows of the matrix U and B.                           //
//     int ncols                                                              //
//        The number of columns of the matrix U.  Also the number of rows and //
//        columns of the matrices D and V.                                    //
//     double* Astar                                                          //
//        On input, a pointer to the first element of an ncols x nrows matrix.//
//        On output, the pseudo-inverse of UDV'.                              //
//                                                                            //
//  Return Values:                                                            //
//        The function is of type void.                                       //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double D[N];                                                           //
//     double Astar[N][M];                                                    //
//     double tolerance;                                                      //
//                                                                            //
//     (your code to initialize the matrices U,D,V)                           //
//                                                                            //
//     Singular_Value_Decomposition_Inverse((double*) U, D, (double*) V,      //
//                                        tolerance, M, N, (double*) Astar);  //
//                                                                            //
//     printf(" The pseudo-inverse of A = UDV' is \n");                       //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //

void Singular_Value_Decomposition_Inverse(double* U, double* D, double* V,
                        double tolerance, int nrows, int ncols, double *Astar)
{
   int i,j,k;
   double *pu, *pv, *pa;
   double dum;

   dum = DBL_EPSILON * D[0] * (double) ncols;
   if (tolerance < dum) tolerance = dum;
   for ( i = 0, pv = V, pa = Astar; i < ncols; i++, pv += ncols)
      for ( j = 0, pu = U; j < nrows; j++, pa++)
        for (k = 0, *pa = 0.0; k < ncols; k++, pu++)
           if (D[k] > tolerance) *pa += *(pv + k) * *pu / D[k];
}
